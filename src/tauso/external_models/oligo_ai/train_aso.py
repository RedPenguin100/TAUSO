import torch
import torch.nn.functional as F
import torch.nn as nn
from torch.optim import Adam
from torch.optim.lr_scheduler import LinearLR

import lightning.pytorch as pl

from lightning.pytorch.loggers.wandb import WandbLogger
from lightning.pytorch.callbacks.lr_monitor import LearningRateMonitor
from lightning.pytorch.callbacks import ModelCheckpoint
from lightning.pytorch.strategies import DDPStrategy

from torchmetrics.regression import R2Score, MeanAbsoluteError, MeanSquaredError

import argparse
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.stats import spearmanr

from src.tauso.external_models.oligo_ai.rinalmo.config import model_config
from src.tauso.external_models.oligo_ai.rinalmo.model.model import RiNALMo
from src.tauso.external_models.oligo_ai.rinalmo.data.alphabet import Alphabet
from src.tauso.external_models.oligo_ai.rinalmo.data.downstream.aso.datamodule import ASODataModule
from src.tauso.external_models.oligo_ai.rinalmo.utils.scaler import StandardScaler
from src.tauso.external_models.oligo_ai.rinalmo.utils.finetune_callback import GradualUnfreezing

class GlobalPooling(nn.Module):
    """Simple and stable global pooling - inspired by RibosomeLoadingPredictionHead"""
    def __init__(self, embed_dim: int, projection_dim: int = 64):
        super().__init__()
        # Optional projection to smaller dimension (can help with computational efficiency)
        self.projection = nn.Linear(embed_dim, projection_dim) if projection_dim < embed_dim else nn.Identity()
        self.output_dim = projection_dim if projection_dim < embed_dim else embed_dim
        
    def forward(self, representations: torch.Tensor, pad_mask: torch.Tensor):
        # Project if needed
        proj_repr = self.projection(representations)  # [batch, seq_len, projection_dim]
        
        # Mask out padded positions
        proj_repr = proj_repr.masked_fill(pad_mask.unsqueeze(-1), 0.0)
        
        # Simple mean pooling over non-padded tokens
        pooled_repr = proj_repr.sum(dim=1) / (~pad_mask).sum(dim=1)[:, None]
        
        return pooled_repr

class ASOPredictionHead(nn.Module):
    """
    Combines ASO and context representations via pooling and concatenation,
    then predicts inhibition using an MLP. Uses method-scaled dosage instead
    of separate dosage and method inputs.
    """
    def __init__(
        self,
        c_in_aso: int,
        c_in_context: int,
        transfection_method_embed_dim: int,
        hidden_dim: int = 128,
        num_layers: int = 3,
        dropout: float = 0.1,
        pooling_projection_dim: int = 64,
    ):
        super().__init__()

        # Use stable global pooling for both ASO and context
        self.aso_pooler = GlobalPooling(c_in_aso, pooling_projection_dim)
        self.context_pooler = GlobalPooling(c_in_context, pooling_projection_dim)

        # The input dimension for the MLP: pooled representations + method-scaled dosage
        mlp_in_dim = self.aso_pooler.output_dim + self.context_pooler.output_dim + transfection_method_embed_dim

        mlp_layers = []
        for i in range(num_layers):
            in_dim = mlp_in_dim if i == 0 else hidden_dim
            out_dim = 1 if i == num_layers - 1 else hidden_dim
            mlp_layers.append(nn.Linear(in_dim, out_dim))
            if i < num_layers - 1:
                mlp_layers.append(nn.ReLU())
                mlp_layers.append(nn.Dropout(dropout))
        self.mlp = nn.Sequential(*mlp_layers)

    def forward(self, aso_repr, context_repr, aso_pad_mask, context_pad_mask, method_scaled_dosage):
        # Pool the ASO and context representations using stable global pooling
        pooled_aso_repr = self.aso_pooler(aso_repr, aso_pad_mask)
        pooled_context_repr = self.context_pooler(context_repr, context_pad_mask)

        # Concatenate the pooled representations and method-scaled dosage
        combined_input = torch.cat([pooled_aso_repr, pooled_context_repr, method_scaled_dosage], dim=-1)

        # Pass through the MLP
        out = self.mlp(combined_input)
        return out.squeeze(-1)

class ASOInhibitionPredictionWrapper(pl.LightningModule):
    def __init__(
        self,
        lm_config: str,
        chem_vocab_size: int,
        chem_embed_dim: int,
        backbone_vocab_size: int,
        backbone_embed_dim: int,
        transfection_method_vocab_size: int,
        transfection_method_embed_dim: int,
        hidden_dim: int,
        num_layers: int,
        dropout: float,
        lr: float,
        scaler: StandardScaler,
        pooling_projection_dim: int = 64,
    ):
        super().__init__()
        self.save_hyperparameters()
        self.scaler = scaler

        self.lm = RiNALMo(model_config(lm_config))
        lm_embed_dim = self.lm.config['model']['transformer'].embed_dim

        self.chem_embedder = nn.Embedding(
            num_embeddings=chem_vocab_size, embedding_dim=chem_embed_dim, padding_idx=0
        )
        self.backbone_embedder = nn.Embedding(
            num_embeddings=backbone_vocab_size, embedding_dim=backbone_embed_dim, padding_idx=0
        )

        # Transfection method embedder - no padding needed for scalar categorical
        self.transfection_method_embedder = nn.Embedding(
            num_embeddings=transfection_method_vocab_size, embedding_dim=transfection_method_embed_dim
        )

        # A small network to better integrate sequence, chemistry, and backbone features
        c_in_aso_combined = lm_embed_dim + chem_embed_dim + backbone_embed_dim
        bottleneck_dim = 128  

        self.aso_feature_combiner = nn.Sequential(
            nn.Linear(c_in_aso_combined, bottleneck_dim),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(bottleneck_dim, lm_embed_dim),
            nn.ReLU()
        )

        self.pred_head = ASOPredictionHead(
            c_in_aso=lm_embed_dim,
            c_in_context=lm_embed_dim,
            transfection_method_embed_dim=transfection_method_embed_dim,
            hidden_dim=hidden_dim,
            num_layers=num_layers,
            dropout=dropout,
            pooling_projection_dim=pooling_projection_dim,
        )

        self.loss = nn.MSELoss()
        self.r2_metric = R2Score()
        self.mae_metric = MeanAbsoluteError()
        self.rmse_metric = MeanSquaredError(squared=False)

        self.lr = lr
        self.pad_idx = self.lm.config['model']['embedding'].padding_idx
        self.val_data_by_custom_id = defaultdict(lambda: {'preds': [], 'targets': []})
        self.test_data_by_custom_id = defaultdict(lambda: {'preds': [], 'targets': []})

    def load_pretrained_lm_weights(self, pretrained_weights_path):
        self.lm.load_state_dict(torch.load(pretrained_weights_path))

    def forward(self, aso_tokens, chem_tokens, backbone_tokens, context_tokens, dosage, transfection_method_tokens):
        aso_repr = self.lm(aso_tokens)["representation"]
        context_repr_original = self.lm(context_tokens)["representation"]

        chem_embeds = self.chem_embedder(chem_tokens)
        backbone_embeds = self.backbone_embedder(backbone_tokens)

        # Scale dosage using log1p before using it
        scaled_dosage = torch.log1p(dosage)

        # Get transfection method embedding and scale by log-transformed dosage
        transfection_method_embeds = self.transfection_method_embedder(transfection_method_tokens)
        # Scale the method embedding by dosage: [batch_size, embed_dim] * [batch_size, 1] -> [batch_size, embed_dim]
        method_scaled_dosage = transfection_method_embeds * scaled_dosage.unsqueeze(-1)

        # Concatenate raw features
        combined_aso_repr_raw = torch.cat([aso_repr, chem_embeds, backbone_embeds], dim=-1)

        # Pass through the feature combiner to get a richer representation
        combined_aso_repr_original = self.aso_feature_combiner(combined_aso_repr_raw)

        aso_pad_mask = aso_tokens.eq(self.pad_idx)
        context_pad_mask = context_tokens.eq(self.pad_idx)

        aso_token_mask = (~aso_pad_mask).unsqueeze(-1).to(combined_aso_repr_original.dtype)
        combined_aso_repr = combined_aso_repr_original * aso_token_mask

        context_token_mask = (~context_pad_mask).unsqueeze(-1).to(context_repr_original.dtype)
        context_repr = context_repr_original * context_token_mask

        # Pass method-scaled dosage to prediction head
        pred = self.pred_head(combined_aso_repr, context_repr, aso_pad_mask, context_pad_mask, method_scaled_dosage)
        return pred

    def _common_step(self, batch, batch_idx, log_prefix: str):
        aso_tokens, chem_tokens, backbone_tokens, context_tokens, inhibition_target, dosage, transfection_method_tokens, custom_ids = batch

        # Scale targets for stable training
        scaled_target = self.scaler.transform(inhibition_target)

        preds_scaled = self(aso_tokens, chem_tokens, backbone_tokens, context_tokens, dosage, transfection_method_tokens)
        loss = self.loss(preds_scaled, scaled_target)

        # For logging and metrics, use the original scale
        preds_unscaled = self.scaler.inverse_transform(preds_scaled.detach())

        self.r2_metric.update(preds_unscaled, inhibition_target)
        self.mae_metric.update(preds_unscaled, inhibition_target)
        self.rmse_metric.update(preds_unscaled, inhibition_target)

        log = { f'{log_prefix}/loss': loss }
        self.log_dict(log, sync_dist=True)
        return loss, preds_unscaled, inhibition_target, custom_ids

    def _eval_step(self, batch, batch_idx, log_prefix, data_collector):
        loss, preds_percent, inhibition_percent, custom_ids = self._common_step(batch, batch_idx, log_prefix=log_prefix)
        for pred, target, custom_id in zip(preds_percent.cpu().numpy(), inhibition_percent.cpu().numpy(), custom_ids):
            data_collector[custom_id]['preds'].append(pred)
            data_collector[custom_id]['targets'].append(target)
        return loss

    def _on_eval_epoch_start(self):
        self.r2_metric.reset(); self.mae_metric.reset(); self.rmse_metric.reset()

    def _calculate_mean_spearman_correlation(self, data_by_custom_id):
        correlations = []
        for custom_id, data in data_by_custom_id.items():
            preds, targets = np.array(data['preds']), np.array(data['targets'])
            if len(preds) > 1:
                corr, _ = spearmanr(preds, targets)
                if not np.isnan(corr): correlations.append(corr)
        return np.mean(correlations) if correlations else 0.0, len(correlations)

    def _calculate_top_pred_target_ratio(self, data_by_custom_id):
        ratios = []
        for custom_id, data in data_by_custom_id.items():
            preds, targets = np.array(data['preds']), np.array(data['targets'])
            if len(preds) > 0:
                top_pred_idx = np.argmax(preds)
                top_pred_target = targets[top_pred_idx]
                mean_target = np.mean(targets)
                if mean_target != 0:
                    ratio = top_pred_target / mean_target
                    ratios.append(ratio)
        return np.median(ratios) if ratios else 0.0, len(ratios)

    def _on_eval_epoch_end(self, log_prefix: str, data_collector):
        if not self.trainer.sanity_checking:
            self.log(f"{log_prefix}/r2", self.r2_metric.compute(), sync_dist=True)
            self.log(f"{log_prefix}/mae_final", self.mae_metric.compute(), sync_dist=True)
            self.log(f"{log_prefix}/rmse", self.rmse_metric.compute(), sync_dist=True)
            mean_spearman, n_custom_ids = self._calculate_mean_spearman_correlation(data_collector)
            self.log(f"{log_prefix}/mean_spearman_corr", mean_spearman, sync_dist=True)
            top_pred_ratio, n_ratios = self._calculate_top_pred_target_ratio(data_collector)
            self.log(f"{log_prefix}/top_pred_target_ratio", top_pred_ratio, sync_dist=True)
            print(f"{log_prefix} - Mean Spearman correlation: {mean_spearman:.4f} across {n_custom_ids} custom_ids")
            print(f"{log_prefix} - Top pred target ratio (median): {top_pred_ratio:.4f} across {n_ratios} custom_ids")
            self.r2_metric.reset(); self.mae_metric.reset(); self.rmse_metric.reset()
        data_collector.clear()

    def training_step(self, batch, batch_idx):
        loss, _, _, _ = self._common_step(batch, batch_idx, log_prefix="train")
        return loss

    def validation_step(self, batch, batch_idx):
        return self._eval_step(batch, batch_idx, log_prefix="val", data_collector=self.val_data_by_custom_id)

    def on_validation_epoch_start(self): return self._on_eval_epoch_start()
    def on_validation_epoch_end(self): return self._on_eval_epoch_end("val", self.val_data_by_custom_id)

    def test_step(self, batch, batch_idx):
        return self._eval_step(batch, batch_idx, log_prefix="test", data_collector=self.test_data_by_custom_id)

    def on_test_epoch_start(self): return self._on_eval_epoch_start()
    def on_test_epoch_end(self): return self._on_eval_epoch_end("test", self.test_data_by_custom_id)

    def configure_optimizers(self):
        optimizer = Adam(filter(lambda p: p.requires_grad, self.parameters()), lr=self.lr)
    
        # Correctly calculate total_steps for the entire training run
        if self.trainer.max_steps and self.trainer.max_steps > 0:
            total_steps = self.trainer.max_steps
        elif self.trainer.max_epochs and self.trainer.max_epochs > 0:
            # Use the dataloader to get the number of batches
            train_dataloader = self.trainer.datamodule.train_dataloader()
            steps_per_epoch = len(train_dataloader)
            total_steps = steps_per_epoch * self.trainer.max_epochs
        else:
            total_steps = self.trainer.estimated_stepping_batches
    
        print(f"Configuring LinearLR scheduler with total_iters = {total_steps}")
        # Consider adding a warm-up phase, which is common for transformers
        # For simplicity, we'll stick to LinearLR but you could change it.
        scheduler = LinearLR(optimizer, start_factor=1.0, end_factor=0.1, total_iters=total_steps)
    
        return {
            "optimizer": optimizer,
            "lr_scheduler": {
                "scheduler": scheduler,
                "interval": "step"
            }
        }

def main(args):
    if args.seed:
        pl.seed_everything(args.seed)
    if args.output_dir:
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    alphabet = Alphabet()
    datamodule = ASODataModule(
        data_path=args.data_path, alphabet=alphabet, batch_size=args.batch_size,
        num_workers=args.num_workers, pin_memory=args.pin_memory, train_ratio=args.train_ratio,
        val_ratio=args.val_ratio, random_state=args.seed if args.seed else 42,
    )
    datamodule.setup(stage='fit')

    # --- Correct Scaler Fitting ---
    # 1. Instantiate the scaler
    scaler = StandardScaler()
    # 2. Extract all target values from the training dataset
    print("Fitting StandardScaler on the training data...")
    train_targets = [item[4] for item in datamodule.train_dataset]
    # 3. Fit the scaler on the entire training set
    scaler.partial_fit(torch.tensor(train_targets, dtype=torch.float32).unsqueeze(-1))

    train_dataset = datamodule.train_dataset.dataset
    chem_vocab_size = len(train_dataset.chem_vocab)
    backbone_vocab_size = len(train_dataset.backbone_vocab)
    transfection_method_vocab_size = len(train_dataset.transfection_method_vocab)

    print(f"Chemistry vocab size: {chem_vocab_size}, Backbone vocab size: {backbone_vocab_size}")
    print(f"Transfection method vocab size: {transfection_method_vocab_size}")

    model = ASOInhibitionPredictionWrapper(
        lm_config=args.lm_config,
        chem_vocab_size=chem_vocab_size,
        chem_embed_dim=args.chem_embed_dim,
        backbone_vocab_size=backbone_vocab_size,
        backbone_embed_dim=args.backbone_embed_dim,
        transfection_method_vocab_size=transfection_method_vocab_size,
        transfection_method_embed_dim=args.transfection_method_embed_dim,
        hidden_dim=args.hidden_dim,
        num_layers=args.num_layers,
        dropout=args.dropout,
        lr=args.lr,
        scaler=scaler,
        pooling_projection_dim=args.pooling_projection_dim,
    )

    if args.pretrained_rinalmo_weights:
        model.load_pretrained_lm_weights(args.pretrained_rinalmo_weights)
    if args.init_params:
        model.load_state_dict(torch.load(args.init_params))

    callbacks = []
    loggers = []
    if args.wandb:
        wandb_logger = WandbLogger(name=args.wandb_experiment_name, save_dir=args.output_dir, project=args.wandb_project, entity=args.wandb_entity, save_code=True)
        loggers.append(wandb_logger)
    if args.checkpoint_every_epoch:
        # Save all epochs (includes best and last)
        callbacks.append(ModelCheckpoint(
            dirpath=args.output_dir,
            filename='aso-epoch_ckpt-{epoch}-{step}',
            every_n_epochs=1,
            save_top_k=-1
        ))
    else:
        # Only save best and last when not saving everything
        callbacks.append(ModelCheckpoint(
            dirpath=args.output_dir,
            filename='aso-best-{epoch}-{val/mae_final:.2f}',
            monitor='val/mae_final',
            mode='min',
            save_top_k=1
        ))
        callbacks.append(ModelCheckpoint(save_last=True))

    if loggers: callbacks.append(LearningRateMonitor(logging_interval="step"))
    if args.ft_schedule: callbacks.append(GradualUnfreezing(unfreeze_schedule_path=args.ft_schedule))

    strategy = 'auto'
    if args.devices != 'auto' and (isinstance(args.devices, str) and "," in args.devices or (isinstance(args.devices, int) and args.devices > 1)):
        strategy = DDPStrategy(find_unused_parameters=False)

    trainer = pl.Trainer(
        accelerator=args.accelerator, devices=args.devices, max_steps=args.max_steps, max_epochs=args.max_epochs,
        gradient_clip_val=args.gradient_clip_val, precision=args.precision, default_root_dir=args.output_dir,
        log_every_n_steps=args.log_every_n_steps, strategy=strategy, logger=loggers, callbacks=callbacks,
    )

    if not args.test_only:
        trainer.fit(model=model, datamodule=datamodule)
    trainer.test(model=model, datamodule=datamodule)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data_path", type=str, help="Path to the CSV file")
    parser.add_argument("--init_params", type=str, default=None, help="Path to model weights")
    parser.add_argument("--output_dir", type=str, default=None, help="Output directory")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument("--checkpoint_every_epoch", action="store_true", default=False)
    parser.add_argument("--test_only", action="store_true", default=False)
    parser.add_argument("--lm_config", type=str, default="giga")
    parser.add_argument("--pretrained_rinalmo_weights", type=str, default=None)
    parser.add_argument("--chem_embed_dim", type=int, default=16, help="Dimension for chemistry embeddings")
    parser.add_argument("--backbone_embed_dim", type=int, default=8, help="Dimension for backbone embeddings")
    parser.add_argument("--transfection_method_embed_dim", type=int, default=4, help="Dimension for transfection method embeddings")
    parser.add_argument("--hidden_dim", type=int, default=128)
    parser.add_argument("--num_layers", type=int, default=3)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--pooling_projection_dim", type=int, default=64, help="Dimension for pooling projection")
    parser.add_argument("--train_ratio", type=float, default=0.8)
    parser.add_argument("--val_ratio", type=float, default=0.1)
    parser.add_argument("--wandb", action="store_true", default=False)
    parser.add_argument("--wandb_experiment_name", type=str, default=None)
    parser.add_argument("--wandb_project", type=str, default=None)
    parser.add_argument("--wandb_entity", type=str, default=None)
    parser.add_argument("--log_every_n_steps", type=int, default=50)
    parser.add_argument("--batch_size", type=int, default=32)
    parser.add_argument("--num_workers", type=int, default=0)
    parser.add_argument("--pin_memory", action="store_true", default=False)
    parser.add_argument("--ft_schedule", type=str, default=None)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--accelerator", type=str, default='auto')
    parser.add_argument("--devices", type=str, default='auto')
    parser.add_argument("--max_steps", type=int, default=-1)
    parser.add_argument("--max_epochs", type=int, default=-1)
    parser.add_argument("--gradient_clip_val", type=float, default=None)
    parser.add_argument("--precision", type=str, default='16-mixed')

    args = parser.parse_args()
    main(args)

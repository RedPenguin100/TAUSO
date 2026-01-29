import torch
import pandas as pd
import numpy as np
from pathlib import Path
from torch.utils.data import DataLoader
import argparse

# TODO: fix
# from .external_models.oligo_ai.rinalmo.data.alphabet import Alphabet
# from .external_models.oligo_ai.rinalmo.utils.scaler import StandardScaler
# from .external_models.oligo_ai.rinalmo.data.downstream.aso.dataset import ASODataset
# from .external_models.oligo_ai.train_aso import ASOInhibitionPredictionWrapper


class ASODatasetWithLen(ASODataset):
    """Wrapper around ASODataset that adds __len__ method for DataLoader compatibility"""

    def __len__(self):
        return len(self.df)


def run_inference(
    model_checkpoint_path: str,
    data_path: str,
    output_path: str = None,
    batch_size: int = 32,
    num_workers: int = 0,
    device: str = "auto"
):
    """
    Run inference on all rows in the dataset and save predictions to CSV.

    Args:
        model_checkpoint_path: Path to the trained model checkpoint
        data_path: Path to the CSV file with splits
        output_path: Path to save the CSV with predictions (optional)
        batch_size: Batch size for inference
        num_workers: Number of workers for data loading
        device: Device to run inference on ('auto', 'cpu', 'cuda')
    """

    # Set device
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")

    # Load the dataframe
    df = pd.read_csv(data_path)

    # Print split statistics
    print("\nDataset split distribution:")
    print(df['split'].value_counts())

    # Load the trained model with its saved scaler
    print(f"Loading model from: {model_checkpoint_path}")
    model = ASOInhibitionPredictionWrapper.load_from_checkpoint(
        model_checkpoint_path
        # Remove the scaler parameter - use the one saved with the model
    )
    model.eval()
    model.to(device)

    print(f"Model loaded successfully with saved scaler: {type(model.scaler)}")

    # Initialize alphabet and dataset
    alphabet = Alphabet()
    dataset = ASODatasetWithLen(
        data_path=data_path,
        alphabet=alphabet,
        pad_to_max_len=True,
    )

    print(f"Loaded dataset with {len(dataset)} samples")

    # Create dataloader for inference
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        num_workers=num_workers,
        pin_memory=True if device == "cuda" else False,
        shuffle=False
    )

    # Run inference
    all_predictions = []

    print("\nStarting inference...")
    with torch.no_grad():
        for batch_idx, batch in enumerate(dataloader):
            aso_tokens, chem_tokens, backbone_tokens, context_tokens, _, dosage, transfection_method_tokens, _ = batch

            # Move tensors to device
            aso_tokens = aso_tokens.to(device)
            chem_tokens = chem_tokens.to(device)
            backbone_tokens = backbone_tokens.to(device)
            context_tokens = context_tokens.to(device)
            dosage = dosage.to(device)
            transfection_method_tokens = transfection_method_tokens.to(device)

            # Use autocast for mixed precision on GPU
            if device == "cuda":
                with torch.autocast(device_type='cuda', dtype=torch.float16):
                    scaled_predictions = model(aso_tokens, chem_tokens, backbone_tokens, context_tokens, dosage, transfection_method_tokens)
            else:
                scaled_predictions = model(aso_tokens, chem_tokens, backbone_tokens, context_tokens, dosage, transfection_method_tokens)

            # Inverse transform the predictions to get the real values
            unscaled_predictions = model.scaler.inverse_transform(scaled_predictions)

            all_predictions.extend(unscaled_predictions.cpu().numpy())

            if (batch_idx + 1) % 10 == 0:
                print(f"Processed {(batch_idx + 1) * batch_size} samples...")

    print(f"Inference completed. Generated {len(all_predictions)} predictions")

    # Verify prediction count matches dataframe
    assert len(all_predictions) == len(df), f"Mismatch: {len(all_predictions)} predictions vs {len(df)} rows"

    # Add predictions to dataframe
    df['predicted_inhibition_percent'] = all_predictions
    df['prediction_error'] = np.abs(np.array(all_predictions) - df['inhibition_percent'].values)
    df['prediction_squared_error'] = (np.array(all_predictions) - df['inhibition_percent'].values) ** 2

    # Calculate overall statistics
    mae = np.mean(df['prediction_error'])
    rmse = np.sqrt(np.mean(df['prediction_squared_error']))

    print(f"\n{'='*60}")
    print(f"Overall Prediction Statistics:")
    print(f"MAE: {mae:.4f}")
    print(f"RMSE: {rmse:.4f}")

    # Calculate R² score
    from scipy.stats import pearsonr, spearmanr
    r, p_value = pearsonr(df['predicted_inhibition_percent'], df['inhibition_percent'])
    r2 = r ** 2
    spearman_corr, _ = spearmanr(df['predicted_inhibition_percent'], df['inhibition_percent'])
    print(f"R²: {r2:.4f}")
    print(f"Pearson correlation: {r:.4f}")
    print(f"Spearman correlation: {spearman_corr:.4f}")

    # Calculate per-split statistics
    print(f"\n{'='*60}")
    print("Per-split Statistics:")
    print(f"{'-'*60}")
    for split in ['train', 'val', 'test']:
        split_df = df[df['split'] == split]
        if not split_df.empty:
            mae = np.mean(split_df['prediction_error'])
            rmse = np.sqrt(np.mean(split_df['prediction_squared_error']))
            r, _ = pearsonr(split_df['predicted_inhibition_percent'], split_df['inhibition_percent'])
            r2 = r ** 2
            spearman_corr, _ = spearmanr(split_df['predicted_inhibition_percent'], split_df['inhibition_percent'])

            print(f"{split.upper():5} - Samples: {len(split_df):6d}")
            print(f"        MAE: {mae:.4f}, RMSE: {rmse:.4f}")
            print(f"        R²: {r2:.4f}, Pearson: {r:.4f}, Spearman: {spearman_corr:.4f}")
            print(f"{'-'*60}")

    # Calculate mean Spearman correlation per custom_id
    print(f"\n{'='*60}")
    print("Group-level Analysis:")
    correlations = []
    for custom_id in df['custom_id'].unique():
        custom_df = df[df['custom_id'] == custom_id]
        if len(custom_df) > 1:
            corr, _ = spearmanr(custom_df['predicted_inhibition_percent'], custom_df['inhibition_percent'])
            if not np.isnan(corr):
                correlations.append(corr)

    if correlations:
        mean_spearman = np.mean(correlations)
        print(f"Mean Spearman correlation across {len(correlations)} custom_ids: {mean_spearman:.4f}")

    # Calculate top prediction target ratio
    ratios = []
    for custom_id in df['custom_id'].unique():
        custom_df = df[df['custom_id'] == custom_id]
        if len(custom_df) > 0:
            top_pred_idx = custom_df['predicted_inhibition_percent'].idxmax()
            top_pred_target = custom_df.loc[top_pred_idx, 'inhibition_percent']
            mean_target = custom_df['inhibition_percent'].mean()
            if mean_target != 0:
                ratio = top_pred_target / mean_target
                ratios.append(ratio)

    if ratios:
        median_ratio = np.median(ratios)
        print(f"Top pred target ratio (median): {median_ratio:.4f} across {len(ratios)} custom_ids")

    # Determine output path
    if output_path is None:
        data_path_obj = Path(data_path)
        output_path = data_path_obj.parent / f"{data_path_obj.stem}.with_predictions.csv"

    # Save results
    df.to_csv(output_path, index=False)
    print(f"\n{'='*60}")
    print(f"Predictions saved to: {output_path}")

    return df


def main():
    parser = argparse.ArgumentParser(description="Run inference on ASO inhibition prediction model")
    parser.add_argument("data_path", type=str, help="Path to the CSV file (preferably with splits)")
    parser.add_argument("--model_checkpoint", type=str, default="../mae_final=19.03.ckpt",
                        help="Path to the trained model checkpoint")
    parser.add_argument("--output_path", type=str, default=None,
                        help="Path to save the CSV with predictions (default: input_file.with_predictions.csv)")
    parser.add_argument("--batch_size", type=int, default=32,
                        help="Batch size for inference")
    parser.add_argument("--num_workers", type=int, default=0,
                        help="Number of workers for data loading")
    parser.add_argument("--device", type=str, default="auto",
                        choices=["auto", "cpu", "cuda"],
                        help="Device to run inference on")

    args = parser.parse_args()

    # Check if files exist
    if not Path(args.data_path).exists():
        raise FileNotFoundError(f"Data file not found: {args.data_path}")

    if not Path(args.model_checkpoint).exists():
        raise FileNotFoundError(f"Model checkpoint not found: {args.model_checkpoint}")

    # Run inference
    results_df = run_inference(
        model_checkpoint_path=args.model_checkpoint,
        data_path=args.data_path,
        output_path=args.output_path,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        device=args.device
    )

    print("\nInference completed successfully!")


if __name__ == "__main__":
    main()
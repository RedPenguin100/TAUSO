**Figure 1. Study overview.**
(a) A gapmer antisense oligonucleotide (ASO) hybridizes to its target transcript and recruits RNase H1, which
cleaves the RNA:DNA duplex to knock down the target. Knockdown efficacy is governed by several competing
determinants — hybridization affinity, local target accessibility (RNA secondary structure), off-target
hybridization, RNase-H1 recruitment, cellular delivery/uptake, and sequence-related toxicity liabilities.
TAUSO learns these signals to rank candidate ASOs by predicted knockdown.
(b) Pipeline. Data: 178,624 ASO–target measurements from the ASO Atlas, filtered to 143,904 strict gapmers,
split by gene into train/test; each ASO is described by six feature classes. Model: XGBoost gradient-boosted
trees trained on a within-experiment ranking objective with gene-grouped cross-validation; 485 per-ASO features and
tuned, regularized parameters define the final model. Application: ranking and design of candidate ASOs
(e.g. DIAPH3), benchmarking against OligoAI, and molecular-dynamics refinement of the top candidates.

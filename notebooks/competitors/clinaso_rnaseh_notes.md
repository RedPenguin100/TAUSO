# ClinASO RNase H1 scoring — notes & learning TODOs

Study summary of what ClinASO's potency score actually is, how it compares to
TAUSO's own RNase H1 features, and what to go learn properly.

## What ClinASO computes (its ONLY potency signal)

`ClinASO_score` = `total_pf` from `ClinASO/asodesigner/scr/_4.py`, transcribed
faithfully in `clinaso_utils.py` (verified ≡ their loop to ~1e-17).

Algorithm: for `offset in 0..3`, take the 13-nt window `ASO[2+offset : 2+offset+13]`,
reverse-complement it (→ the target RNA the ASO hybridizes), look up each base's
weight at its column in a 4×13 **mononucleotide** PWM
(`ClinASO/asodesigner/text/human_rnaseH_pwm.txt`), and sum. Sum the 4 windows → `total_pf`.
Lower is ranked better (sorted ascending in `_8.py`).

## Why it's weak

- It's a plain **linear sum** of PWM weights — **not an occupancy**: no
  exp/Boltzmann weighting, no partition function, unbounded (can't be a binding probability).
- The 4 windows overlap by 10 of 13 nt → the same target base is scored 4× against
  4 *different* PWM columns and summed. Ad hoc smear, not a clean single-site score
  (which would take the max) nor independent sites.
- **No target accessibility / folding ΔG**, no duplex Tm (only crude GC%), no
  chemistry / gapmer-modification modeling, no trained efficacy model.
- The rest of ClinASO (`_5`..`_8`) is off-target / SNP / homology **specificity**
  filtering — that's its real value proposition, not efficacy prediction.

## ClinASO's RNase H1 motifs are NOT the ones we use

| | TAUSO (our study) | ClinASO |
|---|---|---|
| Source | **Kiełpiński et al. 2017**, NAR (PMC5728404), exps R4a/R4b/R7 | undocumented `human_rnaseH_pwm.txt` |
| Encoding | **dinucleotide** (16×W) primary + mononuc + krel/logFC variants | mononucleotide 4×13 only |
| Chemistry | **gap-aware**: scores strictly within the longest DNA gap; modified flanks zeroed | **chemistry-blind**: fixed windows at offsets 2–5 |
| Aggregation | **max** over valid windows (best cut site) | **sum** of 4 overlapping windows |
| Matrix | scissile-site columns hard-zeroed | no zeroing; weight values unrelated to ours |

Conclusion: ClinASO's RNase H1 handling is strictly cruder than TAUSO's.

## Learning TODOs

- [ ] Read **Kiełpiński et al. 2017** (NAR, PMC5728404) — understand the R4a/R4b/R7
      experiments and how the krel (relative cleavage rate) vs logFC weights were derived.
- [ ] Internalize why **dinucleotide context** matters for RNase H1 cleavage
      (the scissile dinucleotide), i.e. why mononucleotide PWMs leave signal on the table.
- [ ] Be able to explain **occupancy** (Boltzmann / partition-function, bounded probability)
      vs a **raw PWM sum**, and why ClinASO's `total_pf` is the latter.
- [ ] Walk through TAUSO's **chemistry-aware, gap-constrained** scoring in
      `src/tauso/features/rnase_motifs/` (`get_longest_dna_gap`, why the 2′-modified
      flanks are zeroed because RNase H1 only cleaves the DNA:RNA duplex).
- [ ] List the axes a real ASO **potency** model needs that ClinASO lacks: target
      accessibility / folding ΔG, ASO:RNA duplex Tm, and modification chemistry.

# Parsimony — held-out test

Within-experiment (custom_id) median Spearman on the held-out test, restricted to the Fig 2 comparison subset (strict 5-10-5 2'-MOE / 3-10-3 cEt gapmers, rows OligoAI also scored). For each K the champion configuration is retrained on train+val using only the top-K features ranked by **gain** importance, then scored once on the test split. Dashed lines are the competitors on the same rows: OligoAI (0.45) and the strongest rule-based tool, OligoWalk·intra (0.17). ★ marks the shipped model, which uses all 485 features (0.58).

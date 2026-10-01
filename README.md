# D-KEFS Analysis

R analysis for the D-KEFS Color-Word Interference Test study (PSY 410 Experiment I).

## Design

One-way repeated-measures ANOVA comparing completion times (seconds) across four conditions:

| Condition | Mean (s) |
|---|---|
| Word Reading | 21.64 |
| Color Naming | 31.27 |
| Inhibition | 52.80 |
| Inhibition/Switching | 59.12 |

## Verified results

- n = 15 complete participants
- F(3, 42) = 46.98, p = 1.77e-13, partial eta-squared = .77
- Mauchly's W = .1685, p = .000393 (sphericity violated)
- Greenhouse-Geisser epsilon = .5255, corrected p < .001
- Bonferroni post-hocs: every pairwise comparison significant except Inhibition vs. Inhibition/Switching (adjusted p = 1.00)

Note: the Greenhouse-Geisser epsilon must be computed from the covariance matrix of
orthonormal contrasts to reproduce .5255.

## Files

- `dkefs_repeated_measures_anova.R` — analysis script (base R only, no extra packages). Reconstructed to reproduce the verified results; run it yourself and confirm the printed values match.
- `Fall_26_Experiment_1_Data_Set_Sheet1.csv` — raw dataset
- `condition-means.png` — condition means plot
- `Ta_PSY410_Experiment1.tex` — LaTeX edition of the course paper (APA student manuscript, `apa7` class). Prose matches the submitted Word manuscript; includes Table 1, Figure 1 (`condition-means.png`), and references. Compile with `pdflatex` (twice, for cross-references); requires the `apa7` document class.

## Usage

Set your working directory to this folder, then `source("dkefs_repeated_measures_anova.R")`.

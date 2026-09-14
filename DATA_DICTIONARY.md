# Aggregate output dictionary

The public CSV files contain one row per drug or one row per prespecified flow stage. They contain no patient, admission, ICU-stay, ECG, or medication-event identifiers.

## phase3b_hc3_precision_results.csv

- `ingredient`: generic drug name.
- `role`: candidate or negative control.
- `n_all`: patients with an informative primary contrast.
- `complete_cases`: patients with all primary-model variables.
- `model_status`: fitted or insufficient complete cases.
- `unadjusted_all_ms`: unadjusted mean among all informative patients.
- `unadjusted_cc_ms`: unadjusted mean in the complete-case sample.
- `hc3_adjusted_ms`, `hc3_se`, `hc3_p`, `hc3_ci_low`, `hc3_ci_high`: primary adjusted estimate and HC3 inference.
- `hc3_fdr_q`: Benjamini-Hochberg q value within the applicable drug-role family.
- `precision_weighted_*`: secondary precision-weighted estimate and inference.
- `weight_p99`, `weight_max`: drug-specific weight diagnostics before capping.
- `hc3_calibrated_*`: secondary empirically calibrated estimate and inference.

## phase3b_replication_check.csv

Drug-level comparison of the original conventional-OLS output with the independently reproduced OLS and HC3 results. `ols_difference_vs_original` should be numerically zero up to floating-point precision.

## drug_flow_counts.csv

Counts and labels used for Supplementary Figure S11. These are aggregate drug-selection counts, not patient counts.

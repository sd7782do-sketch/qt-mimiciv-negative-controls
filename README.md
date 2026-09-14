# Drug-associated QTcF differences in MIMIC-IV

This repository contains code and nonidentifying aggregate outputs supporting the manuscript **Drug-associated differences in the Fridericia-corrected QT interval in critically ill adults: a within-patient analysis of MIMIC-IV with prespecified negative controls**.

## Study overview

The study links MIMIC-IV intensive care medication administrations to MIMIC-IV-ECG measurements. Each patient contributes exposed and unexposed electrocardiograms. The prespecified primary analysis uses Fridericia-corrected QT (QTcF), narrow-QRS electrocardiograms, a 24-hour exposure window, covariate-adjusted ordinary least squares models with HC3 standard errors, Benjamini-Hochberg correction, and eight negative-control drugs. A precision-weighted analysis and empirical calibration are secondary analyses.

## Public contents

- `code/analyse_phase3b.py`: reproduces the primary HC3, precision-weighted, and calibrated estimates from the restricted patient-level analysis input.
- `code/make_release_figures.py`: recreates Figures 2 and 3 and Supplementary Figure S11 from aggregate results and prespecified flow counts.
- `code/validate_release.py`: checks numerical and disclosure safeguards in the public release.
- `results/phase3b_hc3_precision_results.csv`: drug-level aggregate estimates.
- `results/phase3b_replication_check.csv`: aggregate comparison with the original conventional-OLS results.
- `results/drug_flow_counts.csv`: aggregate drug-selection flow counts.
- `figures/`: publication-ready figures and the supplementary flow diagram.
- `docs/`: analysis protocol summary, aggregate-output dictionary, transparency statement, and STROBE/RECORD crosswalk.

## Restricted source data

MIMIC-IV and MIMIC-IV-ECG are credentialed-access resources distributed through PhysioNet. Patient-level, stay-level, ECG-level, and medication-administration-level derived files are deliberately excluded from this release. Eligible researchers must complete the required training, sign the applicable data use agreement, obtain both source datasets from PhysioNet, and reconstruct the analytic input locally.

The Phase 3B verification script expects the restricted input file at:

```text
data/phase3b_patient_covariate_pairs.csv
```

It also expects the aggregate conventional-OLS output at:

```text
data/phase3b_primary_results.csv
```

Do not commit either restricted input file to a public repository.

## Software

The source pipeline used R 4.5.0, data.table 1.18.4, and DuckDB 1.5.5. The independent HC3 verification script uses Python 3 with NumPy, pandas, SciPy, and Matplotlib; exact compatible dependencies are listed in `environment.yml`.

## Reproduction

After placing the authorized restricted inputs under `data/`:

```bash
python code/analyse_phase3b.py
python code/make_release_figures.py
python code/validate_release.py
```

The scripts write aggregate outputs only. Before any public upload, run the validation script and confirm that no file contains patient or stay identifiers.

## Citation

Use the metadata in `CITATION.cff`. Add the GitHub release URL and Figshare DOI only after those records have been created. No DOI is assigned in this package.

## License and data-use restrictions

Repository code is released under the MIT License. This license does not apply to MIMIC-IV, MIMIC-IV-ECG, or any patient-level derivative. Access to and use of those data remain governed by PhysioNet's credentialed data license and data use agreement.

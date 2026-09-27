# QT within-patient reanalysis (MIMIC-IV / MIMIC-IV-ECG)

Code and non-identifying aggregate outputs for a within-patient analysis of the
Fridericia-corrected QT interval in intensive care, including a time-based
pseudo-exposure contrast used to characterise the residual displacement detected
by prespecified negative-control drugs.

No patient-level or stay-level MIMIC-derived data are included. MIMIC-IV and
MIMIC-IV-ECG require PhysioNet credentialing and acceptance of the applicable
data use agreement.

## Contents

### `code/`

| Script | Purpose |
| --- | --- |
| `analyse_phase3b.py` | Original feasibility and primary analysis pipeline (phases 2–3b) |
| `qt_01_extract_admin.R` | Extraction of administration timestamps from eMAR and ICU infusion records (DuckDB). Drug-name patterns are treated as regular expressions, as in the prespecified drug map |
| `qt_02_reanalysis.R` | Exposure classification under look-back, washout and pre-exposure definitions; contrast construction; marginal and HC3 inference; time-based pseudo-exposure contrast; baseline-QTcF stratification; empirical calibration; false discovery rate |
| `qt_03_finish.R` | Coverage diagnostics for the event-based exposure definition; baseline-slope figure; calibration tables; word counts; assembly of this repository |
| `qt_04_fix_missing_drugs.R` | Diagnostic for drugs yielding no extracted administration events |
| `qt_05_acetaminophen_sensitivity.R` | Acetaminophen restricted to single-ingredient products, and refitting of the systematic-error distribution |
| `qt_06_lab_gap.R` | Distribution of the interval between the laboratory measurement and the electrocardiogram |

### `results/`

Aggregate outputs, one file per analysis block. Files numbered `01`–`08` are
produced by `qt_02_reanalysis.R`; files prefixed `A`–`E` by `qt_03`, `qt_05` and
`qt_06`. `feas_00_drug_map.csv` is the prespecified map of candidate drugs,
negative controls and matching patterns. `qt_admin_coverage.csv` and
`qt_emar_event_audit.csv` document the administration extraction.
`sessionInfo.txt` records the computational environment.

### `figures/`

Figures as submitted, at 600 dpi.

### `docs/`

STROBE and RECORD reporting checklists.

## Reproducing

Set the paths in the `CFG` block at the top of each script and run them in
numerical order. `qt_01` must precede `qt_02` if the event-based exposure
definitions are required; the look-back definitions run without it. Requires R
with `data.table`, `ggplot2`, `DBI` and `duckdb`, and credentialed access to
MIMIC-IV and MIMIC-IV-ECG.

## Note on an earlier release

An earlier release of this repository listed directories that were not included
in the archive. This release supersedes it.

## Licence

See `LICENSE`.

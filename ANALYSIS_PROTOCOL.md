# Analysis protocol summary

## Objective

Estimate within-patient associations between documented administration of candidate QT-active drugs and QTcF, while quantifying residual systematic error using prespecified negative controls.

## Design

- Data sources: MIMIC-IV ICU medication administrations linked to MIMIC-IV-ECG.
- Unit of analysis: patient-level exposed-minus-unexposed contrast after averaging informative ICU-stay contrasts within patient.
- Primary outcome: Fridericia-corrected QT interval in milliseconds.
- Primary ECG stratum: QRS duration below 120 ms.
- Primary exposure window: documented administration during the preceding 24 hours.
- Minimum model support: 80 complete patient-drug contrasts.

## Quality filters

- RR interval: 300-2000 ms.
- QT interval: 200-700 ms.
- QRS duration: 40-220 ms.
- QTcF: 250-700 ms.
- Records evaluated: 799,847.
- Records retained after all filters: 797,777.
- Records excluded: 2,070.

## Primary model

Ordinary least squares regression of the patient-level QTcF contrast on exposed-minus-unexposed differences in time from ICU admission, heart rate, potassium, magnesium, total calcium, creatinine, and log(1 + mean ECG count per stay). Covariate differences are mean-centered. With an intercept and mean-centered covariates, the adjusted point estimate equals the complete-case unadjusted mean by construction. HC3 standard errors support confidence intervals and hypothesis tests. Benjamini-Hochberg correction is applied across estimable candidate-drug models.

## Negative controls and calibration

Eight prespecified negative-control drugs characterize residual systematic error. Their HC3 estimates are modeled as normally distributed with fitted mean and between-estimate standard deviation. Candidate-drug estimates are shifted by the fitted mean and their uncertainty incorporates both the individual standard error and fitted systematic-error variance. Calibration is secondary because no synthetic positive controls are available and constant systematic error across true effects is assumed.

## Precision-weighted sensitivity analysis

Weights are proportional to `n_exposed * n_unexposed / (n_exposed + n_unexposed)`, capped at the drug-specific 99th percentile and normalized to mean 1. Covariates are weighted-mean centered and HC3 inference is retained. Repeated ECGs are not treated as independent outcome observations.

## Interpretation

All estimates are associative. The design removes time-invariant between-patient confounding but remains susceptible to time-varying confounding, confounding by indication, informative monitoring, co-medication, reverse causation, exposure-timing error, and complete-case selection.

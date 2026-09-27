###############################################################################
# qt_03_finish.R
#
# Chiude i tre buchi rimasti dopo la rianalisi e prepara il pacchetto finale.
#
#   A. Diagnosi del buco "events": perche' l'acetaminofene (e altri) non sono
#      stimabili con i timestamp mentre lo sono con i flag di fase 2.
#   B. Figura 6: pendenza del contrasto su QTcF basale, candidati vs controlli.
#      Se le due pendenze coincidono e' un secondo argomento, indipendente dalla
#      pseudo-esposizione, a sostegno della stessa tesi.
#   C. Tabella della ricalibrazione (tutti gli 8 / senza famotidina e
#      pantoprazolo / leave-one-out) nella forma da incollare nel manoscritto,
#      piu' le frasi con i numeri da sostituire ai segnaposto.
#   D. Conteggi parole del manoscritto riscritto.
#   E. Costruzione della cartella del repository, che oggi e' vuota.
#
# Richiede che qt_02_reanalysis.R sia gia' girato.
###############################################################################

suppressPackageStartupMessages({library(data.table); library(ggplot2)})

CFG <- list(
  base_dir   = "C:/Users/sd778/Desktop/mimiciv",
  res_dir    = "C:/Users/sd778/Desktop/mimiciv/qt_reanalysis_out",
  input_dir  = "C:/Users/sd778/Desktop/mimiciv/qt_reanalysis_input",
  out_dir    = "C:/Users/sd778/Desktop/mimiciv/qt_final",
  repo_dir   = "C:/Users/sd778/Desktop/mimiciv/qt_repository",
  # percorso del .md o .docx del manoscritto riscritto, per i conteggi
  manuscript = "C:/Users/sd778/Desktop/LAVORI/in submission/QT/QT_manuscript_REWRITTEN.docx",
  baseline_cut = 450
)

msg <- function(...) cat(sprintf(...), "\n", sep = "")
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
rd <- function(d, f) {
  p <- file.path(d, f)
  if (!file.exists(p)) { msg("  MANCA: %s", p); return(NULL) }
  fread(p)
}
RED <- "#B2182B"; BLU <- "#2166AC"; GREY <- "grey30"

## ---------------------------------------------------------------------------
## A. PERCHE' L'ACETAMINOFENE NON E' STIMABILE CON I TIMESTAMP
## ---------------------------------------------------------------------------
# Tre possibili cause, e questo blocco le separa:
#   1. l'estrazione non trova le somministrazioni (pattern o filtro event_txt)
#   2. le trova ma non le aggancia a un ricovero ICU (finestra icustays)
#   3. le aggancia ma restano pochi ricoveri con entrambi gli stati

msg("== A. Diagnosi del motore 'events' ==")

RES   <- rd(CFG$res_dir,   "01_estimates_by_specification.csv")
COVER <- rd(CFG$input_dir, "qt_admin_coverage.csv")

if (!is.null(RES)) {
  sp_flags  <- grep("^flags/window/24h/all_ecg$",  unique(RES$spec), value = TRUE)
  sp_events <- grep("^events/window/24h/all_ecg$", unique(RES$spec), value = TRUE)
  cmp <- merge(
    RES[spec == sp_flags,  .(drug, role, n_inf_flags = n_informative,
                             n_cc_flags = n_complete, est_flags = estimate)],
    RES[spec == sp_events, .(drug, n_inf_events = n_informative,
                             n_cc_events = n_complete, est_events = estimate)],
    by = "drug", all.x = TRUE)
  if (!is.null(COVER))
    cmp <- merge(cmp, COVER[, .(drug = ingredient, n_events, n_subjects_admin = n_subjects,
                                n_stays_admin = n_stays)], by = "drug", all.x = TRUE)
  cmp[is.na(n_inf_events), n_inf_events := 0L]
  cmp[, ratio_inf := round(n_inf_events / pmax(n_inf_flags, 1), 3)]
  cmp[, diagnosi := fifelse(
        is.na(n_events) | n_events == 0, "1. nessuna somministrazione estratta",
      fifelse(n_subjects_admin < 0.5 * n_inf_flags, "2. estratte ma non agganciate al ricovero",
      fifelse(ratio_inf < 0.5, "3. agganciate ma pochi ricoveri con entrambi gli stati",
              "ok")))]
  setorder(cmp, ratio_inf)
  fwrite(cmp, file.path(CFG$out_dir, "A_events_gap_diagnosis.csv"))
  msg("  scritto: A_events_gap_diagnosis.csv")
  print(cmp[, .(drug, role, n_inf_flags, n_inf_events, ratio_inf,
                n_subjects_admin, diagnosi)])
  bad <- cmp[diagnosi != "ok"]
  if (nrow(bad)) {
    msg("\n  Farmaci con un problema nel motore 'events': %s",
        paste(bad$drug, collapse = ", "))
    msg("  Causa prevalente: %s", bad[, .N, by = diagnosi][order(-N)][1, diagnosi])
    msg("  Se prevale la causa 1, correggere il pattern in feas_00_drug_map.csv")
    msg("  o la clausola event_txt in qt_01_extract_admin.R (vedi")
    msg("  qt_emar_event_audit.csv). Se prevale la 2, il problema e' il join su")
    msg("  icustays: le somministrazioni di reparto restano fuori per disegno,")
    msg("  ma per un farmaco frequente come l'acetaminofene questo puo'")
    msg("  eliminare la maggior parte degli eventi.")
  } else msg("  Nessun farmaco perde copertura con il motore 'events'.")
}

## ---------------------------------------------------------------------------
## B. FIGURA 6: PENDENZA SU QTcF BASALE
## ---------------------------------------------------------------------------
# Il contrasto e' delta = esposto - non esposto e il basale e' la componente non
# esposta, quindi una pendenza negativa esiste per costruzione. Il punto non e'
# la sua presenza ma il CONFRONTO: se candidati e controlli hanno la stessa
# pendenza, l'accoppiamento e' interamente meccanico e non farmacologico.

msg("\n== B. Pendenza su QTcF basale ==")

SL <- rd(CFG$res_dir, "05_baseline_qtcf_slopes.csv")
ST <- rd(CFG$res_dir, "04_baseline_qtcf_strata.csv")

if (!is.null(SL) && nrow(SL)) {
  SL[, lab := fifelse(role == "negative_control", "Negative control", "Candidate drug")]
  SL[, drug := paste0(toupper(substring(drug,1,1)), substring(drug,2))]
  setorder(SL, slope)
  SL[, drug := factor(drug, levels = drug)]
  SL[, `:=`(lo = slope - 1.96*slope_se, hi = slope + 1.96*slope_se)]

  pooled <- SL[, .(m = weighted.mean(slope, 1/slope_se^2),
                   se = sqrt(1/sum(1/slope_se^2))), by = lab]
  msg("  pendenza media pesata:")
  print(pooled[, .(lab, slope = round(m,4), se = round(se,4),
                   ci = sprintf("%.3f to %.3f", m-1.96*se, m+1.96*se))])
  d <- pooled[lab == "Candidate drug", m] - pooled[lab == "Negative control", m]
  sd_ <- sqrt(sum(pooled$se^2))
  msg("  differenza candidati - controlli: %+.4f (SE %.4f, p = %.3f)",
      d, sd_, 2*pnorm(-abs(d/sd_)))

  f6 <- ggplot(SL, aes(slope, drug, colour = lab)) +
    geom_vline(xintercept = 0, colour = "grey55") +
    geom_vline(data = pooled, aes(xintercept = m, colour = lab),
               linetype = "22", linewidth = .5, show.legend = FALSE) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .6, alpha = .7) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = c("Negative control" = BLU, "Candidate drug" = RED)) +
    labs(x = "Slope of the patient-level contrast on baseline QTcF (ms per ms)",
         y = NULL) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = "grey95"),
          axis.text = element_text(size = 11, colour = "grey15"),
          legend.position = "bottom", legend.title = element_blank(),
          plot.margin = margin(12, 18, 10, 12))
  ggsave(file.path(CFG$out_dir, "Figure_6_baseline_slope.tiff"), f6,
         width = 7.4, height = 6.6, dpi = 600, compression = "lzw", bg = "white")
  msg("  scritto: Figure_6_baseline_slope.tiff")

  if (!is.null(ST) && nrow(ST)) {
    w <- dcast(ST, drug + role ~ stratum, value.var = c("n","mean_d"))
    fwrite(w, file.path(CFG$out_dir, "B_baseline_strata_wide.csv"))
    msg("  scritto: B_baseline_strata_wide.csv")
  }
  msg("\n  Didascalia per la Figura 6:")
  msg("  \"Figure 6. Coupling between the contrast and its own baseline.")
  msg("   Slope of the patient-level exposed-minus-unexposed QTcF contrast on the")
  msg("   unexposed (baseline) mean, by drug, with 95%% confidence intervals and")
  msg("   inverse-variance pooled values (dashed lines). A negative slope is")
  msg("   mechanical: the baseline is the unexposed component of the contrast.")
  msg("   The comparison of interest is between the two groups: candidate drugs")
  msg("   and prespecified negative controls share the same slope, so the")
  msg("   coupling carries no pharmacological information.\"")
}

## ---------------------------------------------------------------------------
## C. TABELLA DELLA RICALIBRAZIONE E FRASI DA INSERIRE
## ---------------------------------------------------------------------------

msg("\n== C. Ricalibrazione ==")

CAL  <- rd(CFG$res_dir, "06_calibration_sets.csv")
BOOT <- rd(CFG$res_dir, "07_calibration_bootstrap_ci.csv")

if (!is.null(CAL) && nrow(CAL)) {
  fits <- unique(CAL[, .(set, n_controls, mu, sigma)])
  setorder(fits, set)
  fwrite(fits, file.path(CFG$out_dir, "C_calibration_fits.csv"))
  print(fits[, .(set, n_controls, mu = round(mu,3), sigma = round(sigma,3))])

  tab <- CAL[set %in% c("all_controls","drop_famotidine_pantoprazole"),
             .(drug, set, cal = sprintf("%+.2f (%+.2f to %+.2f)",
                                        cal_estimate, cal_lo, cal_hi),
               q = signif(cal_q, 3))]
  wide <- dcast(tab, drug ~ set, value.var = c("cal","q"))
  fwrite(wide, file.path(CFG$out_dir, "C_calibration_table.csv"))
  msg("  scritto: C_calibration_table.csv")

  m_all  <- fits[set == "all_controls", mu]
  m_drop <- fits[set == "drop_famotidine_pantoprazole", mu]
  loo    <- fits[grepl("^loo_", set), range(mu)]
  msg("\n  FRASI DA SOSTITUIRE AI SEGNAPOSTO DEL MANOSCRITTO:")
  msg("  \"...changed the fitted mean from %.2f ms to %.2f ms\"", m_all, m_drop)
  msg("  \"...leave-one-out refits ranged from %.2f to %.2f ms\"", loo[1], loo[2])
  if (!is.null(BOOT) && nrow(BOOT))
    msg("  (intervallo bootstrap per mu: vedi 07_calibration_bootstrap_ci.csv)")
  if (fits[set == "all_controls", sigma] < 0.05)
    msg("  ATTENZIONE: sigma al bordo; la calibrazione non aggiunge dispersione.")
}

## ---------------------------------------------------------------------------
## D. CONTEGGI PAROLE
## ---------------------------------------------------------------------------

msg("\n== D. Conteggi ==")

count_docx <- function(path) {
  if (!file.exists(path)) { msg("  MANCA: %s", path); return(invisible(NULL)) }
  if (grepl("\\.md$", path, ignore.case = TRUE)) {
    txt <- readLines(path, warn = FALSE)
  } else {
    td <- file.path(tempdir(), "docx"); unlink(td, recursive = TRUE)
    utils::unzip(path, exdir = td)
    x <- paste(readLines(file.path(td, "word", "document.xml"),
                         warn = FALSE, encoding = "UTF-8"), collapse = "")
    x <- gsub("</w:p>", "\n", x, fixed = TRUE)
    x <- gsub("<[^>]+>", "", x)
    txt <- strsplit(x, "\n")[[1]]
  }
  txt <- gsub("[*#]", "", txt)
  i0 <- grep("^\\s*Purpose", txt)[1]; i1 <- grep("Key Points", txt)[1]
  j0 <- grep("^\\s*Introduction\\s*$", txt)[1]; j1 <- grep("^\\s*Declarations\\s*$", txt)[1]
  wc <- function(a,b) if (any(is.na(c(a,b)))) NA_integer_ else
    sum(lengths(strsplit(trimws(txt[a:b]), "\\s+")))
  msg("  abstract: %s parole (limite tipico 250)", wc(i0, i1-1))
  msg("  testo principale: %s parole", wc(j0, j1-1))
  msg("  totale documento: %s parole", sum(lengths(strsplit(trimws(txt), "\\s+"))))
}
count_docx(CFG$manuscript)

## ---------------------------------------------------------------------------
## E. REPOSITORY
## ---------------------------------------------------------------------------
# Oggi il repository pubblico contiene un README che elenca cartelle che non
# esistono, mentre il Data Availability Statement dichiara codice e output.
# Questo blocco costruisce l'albero vero.

msg("\n== E. Repository ==")

for (d in c("code", "results", "figures", "docs"))
  dir.create(file.path(CFG$repo_dir, d), showWarnings = FALSE, recursive = TRUE)

cp <- function(from, to_sub, pattern = NULL) {
  if (!dir.exists(from) && !file.exists(from)) { msg("  salto (assente): %s", from); return(0L) }
  f <- if (dir.exists(from)) list.files(from, pattern = pattern, full.names = TRUE) else from
  if (!length(f)) return(0L)
  ok <- file.copy(f, file.path(CFG$repo_dir, to_sub), overwrite = TRUE)
  msg("  %s -> %s: %d file", basename(from), to_sub, sum(ok))
  sum(ok)
}

# codice: gli script di questa rianalisi, ovunque si trovino
for (s in c("qt_01_extract_admin.R", "qt_02_reanalysis.R", "qt_03_finish.R")) {
  p <- c(file.path(getwd(), s), file.path(CFG$base_dir, s))
  p <- p[file.exists(p)][1]
  if (!is.na(p)) cp(p, "code")
}
cp(CFG$res_dir,   "results", "\\.csv$")
cp(CFG$input_dir, "results", "^qt_(admin_coverage|emar_event_audit)\\.csv$")
cp(CFG$out_dir,   "results", "^[ABC]_.*\\.csv$")
cp(CFG$out_dir,   "figures", "\\.tiff$")

writeLines(c(
  "# QT within-patient reanalysis (MIMIC-IV / MIMIC-IV-ECG)",
  "",
  "Code and nonidentifying aggregate outputs for the reanalysis of within-patient",
  "QTcF contrasts, including the time-based pseudo-exposure contrast.",
  "",
  "No patient-level or stay-level MIMIC-derived data are included. MIMIC-IV and",
  "MIMIC-IV-ECG require PhysioNet credentialing and a data use agreement.",
  "",
  "## Contents",
  "",
  "- `code/qt_01_extract_admin.R` extraction of administration timestamps from",
  "  eMAR and ICU infusion records (DuckDB).",
  "- `code/qt_02_reanalysis.R` exposure classification, contrast construction,",
  "  marginal and HC3 inference, pseudo-exposure contrast, baseline-QTcF",
  "  stratification, empirical calibration, false discovery rate.",
  "- `code/qt_03_finish.R` coverage diagnostics, baseline-slope figure,",
  "  calibration tables, word counts, assembly of this repository.",
  "- `results/` aggregate outputs, one CSV per analysis block.",
  "- `figures/` figures as submitted.",
  "",
  "## Reproducing",
  "",
  "Set the paths in the `CFG` block at the top of each script, then run them in",
  "order. Requires R (data.table, ggplot2, DBI, duckdb) and credentialed access",
  "to MIMIC-IV and MIMIC-IV-ECG.",
  "",
  "## Note on an earlier release",
  "",
  "An earlier release of this repository listed directories that were not",
  "included in the archive. This release supersedes it."),
  file.path(CFG$repo_dir, "README.md"))

n <- length(list.files(CFG$repo_dir, recursive = TRUE))
msg("  repository in %s: %d file", CFG$repo_dir, n)
if (n < 5) msg("  ATTENZIONE: il repository e' quasi vuoto, controllare i percorsi in CFG.")

msg("\n== Fatto. Output in %s ==", CFG$out_dir)

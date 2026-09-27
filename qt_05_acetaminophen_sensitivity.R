###############################################################################
# qt_05_acetaminophen_sensitivity.R
#
# Il pattern della mappa, 'acetaminophen|paracetamol', cattura anche le
# associazioni: nell'eMAR di MIMIC ci sono ~11.200 righe di
# hydrocodone-acetaminophen, ~11.000 di oxycodone-acetaminophen e ~20.500 di
# acetaminophen-caff-butalbital. L'acetaminofene e' il controllo negativo con
# il displacement piu' grande (-4.05 ms su 2.670 pazienti), quindi la presenza
# di oppioidi al suo interno e' l'obiezione piu' ovvia.
#
# Questo script ri-estrae l'acetaminofene limitandolo ai prodotti a principio
# attivo singolo e ricalcola: il contrasto, la pseudo-esposizione temporale e
# la distribuzione dell'errore sistematico con gli altri sette controlli.
#
# La sensibilita' e' possibile SOLO con il motore 'events': i flag
# exposed_6h/12h/24h di fase 2 sono gia' calcolati con il pattern largo e non
# sono modificabili a valle. Va dichiarato cosi' nei Metodi.
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(DBI); library(duckdb)
})
set.seed(20260927)

CFG <- list(
  base_dir = "C:/Users/sd778/Desktop/mimiciv",
  f_emar     = "emar.csv.gz",
  f_icustays = "icustays.csv.gz",
  f_phase2   = "qt_feasibility_phase2/phase2_ecg_drug_long.parquet",
  f_phase3   = "qt_feasibility_phase3/phase3_ecg_lab_nearest.parquet",
  res_dir    = "C:/Users/sd778/Desktop/mimiciv/qt_reanalysis_out",
  out_dir    = "C:/Users/sd778/Desktop/mimiciv/qt_final",

  # stretto: deve COMINCIARE con il principio attivo e non essere seguito da
  # un trattino. Cattura 'acetaminophen', 'acetaminophen iv',
  # 'acetaminophen (liquid)'; esclude 'hydrocodone-acetaminophen' (non in
  # posizione iniziale) e 'acetaminophen-caff-butalbital' (trattino dopo).
  pattern_strict = '^(acetaminophen|paracetamol)([ (]|$)',
  pattern_loose  = 'acetaminophen|paracetamol',

  qrs_narrow_max = 120,
  qtcf_range = c(250, 700),
  exposure_window_h = 24,
  washout_h = 24,
  lab_cols_use = c("potassium","magnesium","calcium..total","creatinine"),
  min_complete_cases = 80
)

msg <- function(...) cat(sprintf(...), "\n", sep = "")
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
pth <- function(f) normalizePath(file.path(CFG$base_dir, f), winslash = "/")

## ---------------------------------------------------------------------------
## 1. RI-ESTRAZIONE CON I DUE PATTERN
## ---------------------------------------------------------------------------

msg("== 1. Estrazione acetaminofene, pattern largo vs stretto ==")

con <- dbConnect(duckdb::duckdb()); on.exit(dbDisconnect(con, shutdown = TRUE))
dbExecute(con, "SET memory_limit='8GB'")

extract_one <- function(pattern, label) {
  q <- sprintf("
    WITH e AS (
      SELECT subject_id, hadm_id,
             CAST(charttime AS TIMESTAMP) AS admin_time,
             lower(medication) AS medication
      FROM read_csv_auto('%s', header = true)
      WHERE medication IS NOT NULL AND charttime IS NOT NULL
        AND (event_txt ILIKE '%%administered%%'
             OR event_txt IN ('Started','Restarted','Applied'))
        AND event_txt NOT ILIKE '%%not given%%'
        AND event_txt NOT ILIKE '%%held%%'
        AND event_txt NOT ILIKE '%%stopped%%'
        AND regexp_matches(lower(medication), '%s')
    ),
    i AS (SELECT subject_id, hadm_id, stay_id,
                 CAST(intime AS TIMESTAMP) intime, CAST(outtime AS TIMESTAMP) outtime
          FROM read_csv_auto('%s', header = true))
    SELECT DISTINCT e.subject_id, i.stay_id, e.admin_time
    FROM e JOIN i
      ON i.subject_id = e.subject_id
     AND (e.hadm_id IS NULL OR i.hadm_id = e.hadm_id)
     AND e.admin_time BETWEEN i.intime AND i.outtime",
    pth(CFG$f_emar), gsub("'", "''", pattern), pth(CFG$f_icustays))
  x <- as.data.table(dbGetQuery(con, q))
  # Chiavi come carattere, come nella tabella ECG: gli id estratti da DuckDB
  # arrivano numerici e il rolling join fallirebbe con "Incompatible join types".
  for (k in intersect(c("subject_id","stay_id"), names(x)))
    set(x, j = k, value = as.character(x[[k]]))
  x[, admin_time := as.POSIXct(admin_time, tz = "UTC")]
  msg("  %-8s eventi: %s | pazienti: %s", label,
      format(nrow(x), big.mark = ","), format(uniqueN(x$subject_id), big.mark = ","))
  x
}

# audit delle stringhe escluse, da riportare nel supplemento
audit <- as.data.table(dbGetQuery(con, sprintf("
  SELECT lower(medication) AS medication, count(*) AS n,
         regexp_matches(lower(medication), '%s') AS kept_strict
  FROM read_csv_auto('%s', header = true)
  WHERE regexp_matches(lower(medication), '%s')
  GROUP BY 1,3 ORDER BY n DESC",
  gsub("'","''",CFG$pattern_strict), pth(CFG$f_emar),
  gsub("'","''",CFG$pattern_loose))))
fwrite(audit, file.path(CFG$out_dir, "D_acetaminophen_product_audit.csv"))
msg("\n  prodotti catturati dal pattern largo:")
print(head(audit, 12))
msg("  righe escluse dal pattern stretto: %s su %s (%.1f%%)",
    format(audit[kept_strict == FALSE, sum(n)], big.mark = ","),
    format(audit[, sum(n)], big.mark = ","),
    100 * audit[kept_strict == FALSE, sum(n)] / audit[, sum(n)])

adm_loose  <- extract_one(CFG$pattern_loose,  "largo")
adm_strict <- extract_one(CFG$pattern_strict, "stretto")

## ---------------------------------------------------------------------------
## 2. ECG E COVARIATE
## ---------------------------------------------------------------------------

read_par <- function(p) as.data.table(dbGetQuery(con,
  sprintf("SELECT * FROM read_parquet('%s')", pth(p))))

msg("\n== 2. ECG e covariate ==")
long <- read_par(CFG$f_phase2)
long[, ecg_time := as.POSIXct(ecg_time, tz = "UTC")]
for (k in c("subject_id","stay_id","study_id"))
  set(long, j = k, value = as.character(long[[k]]))

ecg <- unique(long[, .(subject_id, stay_id, study_id, ecg_time,
                       hours_from_icu_admit, heart_rate, qtcf_ms, qrs_ms)],
              by = "study_id")
ecg[, `:=`(qtcf = as.numeric(qtcf_ms), hr = as.numeric(heart_rate),
           time_from_icu_h = as.numeric(hours_from_icu_admit))]
ecg <- ecg[is.finite(qtcf) & qtcf %between% CFG$qtcf_range &
           is.finite(qrs_ms) & qrs_ms < CFG$qrs_narrow_max]

labs <- read_par(CFG$f_phase3)
labs[, study_id := as.character(study_id)]
labs[, label := make.names(tolower(trimws(label)))]
labs <- labs[label %in% CFG$lab_cols_use]
if ("lab_gap_min" %in% names(labs)) setorder(labs, study_id, label, lab_gap_min)
labs <- unique(labs, by = c("study_id","label"))
for (lb in CFG$lab_cols_use) {
  x <- labs[label == lb, .(study_id, v = as.numeric(lab_value))]
  setnames(x, "v", lb); ecg <- x[ecg, on = "study_id"]
}
msg("  ECG narrow-QRS: %s", format(nrow(ecg), big.mark = ","))

VC <- c("qtcf","hr","time_from_icu_h", CFG$lab_cols_use)

## ---------------------------------------------------------------------------
## 3. CONTRASTO E STIMA
## ---------------------------------------------------------------------------

classify <- function(adm, mode) {
  ad <- adm[, .(subject_id, charttime = admin_time)]
  e <- copy(ecg); setnames(e, "ecg_time", "charttime"); e[, row_id := .I]
  setorder(ad, subject_id, charttime); setorder(e, subject_id, charttime)
  j <- ad[e, on = .(subject_id, charttime), roll = TRUE,
          .(row_id = i.row_id, last_adm = x.charttime)]
  e[, last_adm := j[["last_adm"]][match(row_id, j[["row_id"]])]]
  e[, gap_h := as.numeric(difftime(charttime, last_adm, units = "hours"))]
  e[, state := NA_character_]
  e[!is.na(gap_h) & gap_h >= 0 & gap_h <= CFG$exposure_window_h, state := "exposed"]
  if (mode == "window")
    e[is.na(state) & (is.na(gap_h) | gap_h > CFG$exposure_window_h), state := "unexposed"]
  else if (mode == "washout")
    e[is.na(state) & (is.na(gap_h) | gap_h > CFG$washout_h), state := "unexposed"]
  else
    e[is.na(state) & is.na(gap_h), state := "unexposed"]
  setnames(e, "charttime", "ecg_time"); e[!is.na(state)]
}

contrast <- function(es) {
  cnt <- es[, .(n_ecg_state = .N), by = .(subject_id, stay_id, state)]
  agg <- es[, lapply(.SD, mean, na.rm = TRUE),
            by = .(subject_id, stay_id, state), .SDcols = VC]
  agg <- agg[cnt, on = c("subject_id","stay_id","state"), nomatch = NULL]
  keep <- c("subject_id","stay_id", VC, "n_ecg_state")
  ex <- agg[state == "exposed", ..keep]; un <- agg[state == "unexposed", ..keep]
  if (!nrow(ex) || !nrow(un)) return(NULL)
  setnames(ex, c(VC,"n_ecg_state"), paste0(c(VC,"n_ecg_state"), "_exposed"))
  setnames(un, c(VC,"n_ecg_state"), paste0(c(VC,"n_ecg_state"), "_unexposed"))
  w <- ex[un, on = c("subject_id","stay_id"), nomatch = NULL]
  for (v in VC) set(w, j = paste0("d_", v),
                    value = w[[paste0(v,"_exposed")]] - w[[paste0(v,"_unexposed")]])
  set(w, j = "d_log_ecg_count",
      value = log1p(w[["n_ecg_state_exposed"]]) - log1p(w[["n_ecg_state_unexposed"]]))
  dc <- c(paste0("d_", VC), "d_log_ecg_count")
  w[, lapply(.SD, mean, na.rm = TRUE), by = subject_id, .SDcols = dc]
}

estimate <- function(pat) {
  cv <- c("d_time_from_icu_h","d_hr", paste0("d_", CFG$lab_cols_use), "d_log_ecg_count")
  cc <- pat[complete.cases(pat[, c("d_qtcf", cv), with = FALSE])]
  if (nrow(cc) < CFG$min_complete_cases) return(list(n = nrow(cc), est = NA, se = NA))
  list(n = nrow(cc), est = cc[, mean(d_qtcf)], se = cc[, sd(d_qtcf)/sqrt(.N)])
}

msg("\n== 3. Contrasto acetaminofene ==")
OUT <- rbindlist(lapply(c("window","washout","pre_only"), function(md)
  rbindlist(lapply(c("loose","strict"), function(lb) {
    adm <- if (lb == "loose") adm_loose else adm_strict
    es <- classify(adm, md); pat <- contrast(es)
    if (is.null(pat)) return(NULL)
    r <- estimate(pat)
    data.table(mode = md, pattern = lb, n_complete = r$n,
               estimate = r$est, se = r$se,
               lo = r$est - 1.96*r$se, hi = r$est + 1.96*r$se)
  }))))
fwrite(OUT, file.path(CFG$out_dir, "D_acetaminophen_sensitivity.csv"))
print(OUT[, .(mode, pattern, n_complete, estimate = round(estimate,2),
              ci = sprintf("%+.2f to %+.2f", lo, hi))])

## ---------------------------------------------------------------------------
## 4. RICALIBRAZIONE CON L'ACETAMINOFENE STRETTO
## ---------------------------------------------------------------------------

msg("\n== 4. Errore sistematico con l'acetaminofene ristretto ==")

RES <- fread(file.path(CFG$res_dir, "01_estimates_by_specification.csv"))
ctrl <- RES[spec == "events/window/24h/all_ecg" & role == "negative_control" &
            estimable == TRUE, .(drug, estimate, se = se_marginal)]

fit <- function(y, se) {
  nll <- function(p) { v <- exp(p[2])^2 + se^2
    0.5*sum(log(2*pi*v) + (y - p[1])^2/v) }
  o <- optim(c(mean(y), log(max(sd(y), 1e-3))), nll, method = "BFGS")
  c(mu = o$par[1], sigma = exp(o$par[2]))
}

new <- OUT[mode == "window" & pattern == "strict"]
c2 <- copy(ctrl)
if (nrow(new) && is.finite(new$estimate))
  c2[drug == "acetaminophen", `:=`(estimate = new$estimate, se = new$se)]

# as.list() su un vettore nominato non genera le colonne dentro data.table():
# le si costruisce esplicitamente.
mk_fit <- function(label, y, se) {
  f <- fit(y, se)
  data.table(set = label, n_controls = length(y),
             mu = unname(f[["mu"]]), sigma = unname(f[["sigma"]]))
}
FITS <- rbind(
  mk_fit("acetaminophen, loose pattern",  ctrl$estimate, ctrl$se),
  mk_fit("acetaminophen, strict pattern", c2$estimate,   c2$se),
  mk_fit("excluding acetaminophen",
         ctrl[drug != "acetaminophen", estimate],
         ctrl[drug != "acetaminophen", se]))
fwrite(FITS, file.path(CFG$out_dir, "D_calibration_with_strict_acetaminophen.csv"))
print(FITS[, .(set, n_controls, mu = round(mu,3), sigma = round(sigma,3))])

msg("\n  FRASE PER I RISULTATI:")
msg("  \"Restricting acetaminophen to single-ingredient products (%s of %s eMAR",
    format(audit[kept_strict == TRUE, sum(n)], big.mark = ","),
    format(audit[, sum(n)], big.mark = ","))
msg("   rows) changed its contrast from %+.2f to %+.2f ms and the fitted",
    OUT[mode=="window" & pattern=="loose", estimate],
    OUT[mode=="window" & pattern=="strict", estimate])
msg("   systematic-error mean from %.2f to %.2f ms.\"",
    FITS[1, mu], FITS[2, mu])

msg("\n== Fatto. Output in %s ==", CFG$out_dir)

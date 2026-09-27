###############################################################################
# qt_02_reanalysis.R
#
# Rianalisi del contrasto QTcF within-patient, adattata ALLO SCHEMA REALE dei
# file prodotti dalla pipeline originale.
#
# Input
#   qt_feasibility_phase2/phase2_ecg_drug_long.parquet
#       1.723.896 righe = 71.829 ECG x 24 farmaci (16 candidati + 8 controlli)
#       colonne: subject_id, hadm_id, stay_id, study_id, ecg_time,
#                hours_from_icu_admit, rr_interval, heart_rate, qt_ms, qrs_ms,
#                jt_ms, qtcf_ms, qtcb_ms, jtcf_ms, jtcb_ms, wide_qrs,
#                ingredient, role, exposed_6h, exposed_12h, exposed_24h
#   qt_feasibility_phase3/phase3_ecg_lab_nearest.parquet
#       laboratori in formato LUNGO: study_id, label, lab_value, lab_gap_min
#   qt_admin_events.parquet  (opzionale, da qt_01_extract_admin.R)
#       subject_id, hadm_id, stay_id, ingredient, role, admin_time, source
#
# Due motori di esposizione
#   "flags"  usa exposed_6h/12h/24h come nella pipeline originale.
#            Riproduce l'analisi pubblicata. NON permette washout ne' pre-only:
#            dai booleani non si ricava la distanza dall'ultima dose.
#   "events" usa i timestamp estratti da qt_01. Permette washout e pre-only,
#            cioe' la correzione per il carryover (amiodarone).
#
# Blocchi
#   A input e audit      B classificazione esposizione     C contrasti
#   D SE marginale vs HC3 condizionale + bootstrap
#   E pseudo-esposizione temporale (test di regressione verso la media)
#   F stratificazione per QTcF basale, estesa ai controlli negativi
#   G calibrazione: tutti gli 8, senza famotidina+pantoprazolo, LOO, bootstrap
#   H FDR ricalcolato sulle SE marginali
###############################################################################

suppressPackageStartupMessages(library(data.table))
set.seed(20260926)

## ---------------------------------------------------------------------------
## CONFIGURAZIONE
## ---------------------------------------------------------------------------

CFG <- list(
  base_dir = "C:/Users/sd778/Desktop/mimiciv",

  f_phase2 = "qt_feasibility_phase2/phase2_ecg_drug_long.parquet",
  f_phase3 = "qt_feasibility_phase3/phase3_ecg_lab_nearest.parquet",
  # NULL = solo motore "flags"; altrimenti il parquet prodotto da qt_01
  f_admin  = "qt_reanalysis_input/qt_admin_events.parquet",

  out_dir = "C:/Users/sd778/Desktop/mimiciv/qt_reanalysis_out",

  exposure_engine = "both",     # "flags" | "events" | "both"

  qrs_narrow_max = 120,
  qtcf_range     = c(250, 700),

  # covariate di laboratorio da usare. NULL = tutte quelle presenti in
  # phase3. I nomi sono normalizzati con make.names(), quindi
  # "calcium, total" diventa "calcium..total".
  lab_cols_use = c("potassium", "magnesium", "calcium..total", "creatinine"),

  # distanza massima ammessa fra prelievo e ECG. Inf riproduce l'originale;
  # 240 e' un valore ragionevole da riportare come sensibilita'. Nei primi
  # record di phase3 si vedono gap di oltre 500 minuti.
  max_lab_gap_min = Inf,

  exposure_window_h = 24,       # finestra primaria: 6, 12 o 24
  washout_default_h = 24,
  washout_by_drug = c(
    amiodarone = 720, methadone = 120, fluconazole = 72, azithromycin = 72,
    levofloxacin = 48, ciprofloxacin = 48, haloperidol = 24,
    propofol = 12, dexmedetomidine = 12),

  exposure_level = "subject",   # "subject" | "stay"  (carryover fra ricoveri)
  contrast_level = "stay",      # "stay" | "patient"

  min_complete_cases = 80,
  baseline_qtcf_cut  = 450,

  B_mean  = 2000,
  B_calib = 200,
  n_threads = max(1L, parallel::detectCores() - 1L)
)

setDTthreads(CFG$n_threads)
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(sprintf(...), "\n", sep = "")
wr <- function(x, name) {
  if (is.null(x) || !nrow(x)) { msg("  (vuoto, non scritto: %s)", name); return(invisible(NULL)) }
  fwrite(x, file.path(CFG$out_dir, paste0(name, ".csv")))
  msg("  scritto: %s.csv (%d righe)", name, nrow(x))
  invisible(x)
}

# HC3 calcolato internamente: niente sandwich/lmtest, e niente interazione fra
# data.table e i pacchetti di modellazione.
hc3_vcov <- function(m) {
  X <- model.matrix(m); e <- as.vector(residuals(m))
  h <- stats::hat(X, intercept = FALSE); h[h > 1 - 1e-10] <- 1 - 1e-10
  u <- e / (1 - h)
  XtXinv <- chol2inv(qr.R(qr(X)))
  XtXinv %*% crossprod(X * u) %*% XtXinv
}

## ---------------------------------------------------------------------------
## LETTURA (duckdb -> arrow -> fread)
## ---------------------------------------------------------------------------

read_tbl <- function(path, where = NULL) {
  if (!file.exists(path)) stop("File non trovato: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "parquet") {
    if (requireNamespace("duckdb", quietly = TRUE) &&
        requireNamespace("DBI", quietly = TRUE)) {
      con <- DBI::dbConnect(duckdb::duckdb())
      on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
      q <- sprintf("SELECT * FROM read_parquet('%s')",
                   normalizePath(path, winslash = "/"))
      if (!is.null(where)) q <- paste(q, "WHERE", where)
      return(as.data.table(DBI::dbGetQuery(con, q)))
    }
    if (requireNamespace("arrow", quietly = TRUE))
      return(as.data.table(arrow::read_parquet(path)))
    stop("Per leggere i parquet serve duckdb oppure arrow.")
  }
  as.data.table(data.table::fread(path))
}

pth <- function(f) file.path(CFG$base_dir, f)

## ---------------------------------------------------------------------------
## A. INPUT
## ---------------------------------------------------------------------------

msg("R %s | data.table %s", getRversion(), utils::packageVersion("data.table"))
msg("== A. Input ==")

long <- read_tbl(pth(CFG$f_phase2))
need <- c("subject_id","stay_id","study_id","ecg_time","hours_from_icu_admit",
          "heart_rate","qtcf_ms","qrs_ms","ingredient","role",
          "exposed_6h","exposed_12h","exposed_24h")
miss <- setdiff(need, names(long))
if (length(miss)) stop("phase2: colonne mancanti -> ", paste(miss, collapse = ", "))

long[, ecg_time := as.POSIXct(ecg_time, tz = "UTC")]
for (k in c("subject_id","stay_id","study_id"))
  set(long, j = k, value = as.character(long[[k]]))
long[, ingredient := tolower(trimws(ingredient))]
long[, role := tolower(trimws(role))]

msg("  phase2: %s righe | %d ECG | %d farmaci | %d pazienti",
    format(nrow(long), big.mark = ","), uniqueN(long$study_id),
    uniqueN(long$ingredient), uniqueN(long$subject_id))

ALL_DRUGS  <- long[, sort(unique(ingredient))]
CONTROLS   <- long[role == "negative_control", sort(unique(ingredient))]
CANDIDATES <- setdiff(ALL_DRUGS, CONTROLS)
msg("  %d candidati, %d controlli negativi", length(CANDIDATES), length(CONTROLS))

# tabella a livello di ECG (una riga per study_id)
ECG_COLS <- c("subject_id","stay_id","study_id","ecg_time","hours_from_icu_admit",
              "heart_rate","qtcf_ms","qrs_ms")
ecg <- unique(long[, ..ECG_COLS], by = "study_id")
ecg[, `:=`(qtcf = as.numeric(qtcf_ms),
           hr = as.numeric(heart_rate),
           time_from_icu_h = as.numeric(hours_from_icu_admit))]
ecg <- ecg[is.finite(qtcf) & qtcf %between% CFG$qtcf_range]

## --- laboratori: formato lungo -> largo per study_id ---
labs_long <- read_tbl(pth(CFG$f_phase3))
if (!all(c("study_id","label","lab_value") %in% names(labs_long)))
  stop("phase3: attese le colonne study_id, label, lab_value")
labs_long[, study_id := as.character(study_id)]
labs_long[, label := tolower(trimws(label))]
# I nomi finiscono in colonne e poi in as.formula(): un'etichetta come
# "calcium, total" romperebbe la formula (la virgola separa gli argomenti) e
# i contrasti verrebbero costruiti solo dove il caso degenera.
labs_long[, label := make.names(label)]
if (is.finite(CFG$max_lab_gap_min) && "lab_gap_min" %in% names(labs_long)) {
  n0 <- nrow(labs_long)
  labs_long <- labs_long[is.na(lab_gap_min) | lab_gap_min <= CFG$max_lab_gap_min]
  msg("  lab: %d -> %d righe con gap <= %s min", n0, nrow(labs_long),
      CFG$max_lab_gap_min)
}
if ("lab_gap_min" %in% names(labs_long)) {
  gq <- labs_long[, as.list(round(quantile(lab_gap_min, c(.5,.75,.9,.95,.99),
                                           na.rm = TRUE)))]
  msg("  lab_gap_min (mediana/75/90/95/99): %s",
      paste(unlist(gq), collapse = " / "))
}
# un valore per study_id x label (il piu' vicino se il gap e' disponibile)
if ("lab_gap_min" %in% names(labs_long)) setorder(labs_long, study_id, label, lab_gap_min)
labs_long <- unique(labs_long, by = c("study_id","label"))

LAB_COLS <- labs_long[, sort(unique(label))]
msg("  laboratori presenti: %s", paste(LAB_COLS, collapse = ", "))
if (!is.null(CFG$lab_cols_use)) {
  miss <- setdiff(CFG$lab_cols_use, LAB_COLS)
  if (length(miss))
    stop("lab_cols_use: etichette assenti in phase3 -> ",
         paste(miss, collapse = ", "),
         "\n  presenti: ", paste(LAB_COLS, collapse = ", "))
  LAB_COLS <- CFG$lab_cols_use
  labs_long <- labs_long[label %in% LAB_COLS]
  msg("  laboratori usati: %s", paste(LAB_COLS, collapse = ", "))
}

# reshape esplicito, senza dcast
labs_wide <- unique(labs_long[, .(study_id)])
for (lb in LAB_COLS) {
  x <- labs_long[label == lb, .(study_id, v = as.numeric(lab_value))]
  setnames(x, "v", lb)
  labs_wide <- x[labs_wide, on = "study_id"]
}
ecg <- labs_wide[ecg, on = "study_id"]

ecg_narrow <- ecg[is.finite(qrs_ms) & qrs_ms < CFG$qrs_narrow_max]
msg("  ECG narrow-QRS: %d su %d", nrow(ecg_narrow), nrow(ecg))

## --- eventi di somministrazione (opzionale) ---
admin <- NULL
if (!is.null(CFG$f_admin) && file.exists(pth(CFG$f_admin))) {
  admin <- read_tbl(pth(CFG$f_admin))
  admin[, admin_time := as.POSIXct(admin_time, tz = "UTC")]
  for (k in c("subject_id","stay_id")) set(admin, j = k, value = as.character(admin[[k]]))
  admin[, ingredient := tolower(trimws(ingredient))]
  msg("  eventi di somministrazione: %s (%d farmaci)",
      format(nrow(admin), big.mark = ","), uniqueN(admin$ingredient))
} else {
  msg("  eventi di somministrazione: ASSENTI -> solo motore 'flags'")
}

ENGINES <- switch(CFG$exposure_engine,
  flags  = "flags",
  events = "events",
  both   = if (is.null(admin)) "flags" else c("flags","events"),
  stop("exposure_engine non valido"))
if ("events" %in% ENGINES && is.null(admin))
  stop("exposure_engine richiede gli eventi ma qt_admin_events.parquet non c'e'. ",
       "Eseguire prima qt_01_extract_admin.R.")

## ---------------------------------------------------------------------------
## B. CLASSIFICAZIONE DELL'ESPOSIZIONE
## ---------------------------------------------------------------------------

washout_for <- function(d) {
  w <- CFG$washout_by_drug[d]
  if (is.na(w)) CFG$washout_default_h else as.numeric(w)
}

# motore "flags": legge direttamente exposed_6h/12h/24h dalla tabella lunga
classify_flags <- function(the_drug, window_h = CFG$exposure_window_h) {
  fl <- paste0("exposed_", window_h, "h")
  if (!fl %in% names(long)) stop("colonna assente: ", fl)
  x <- long[ingredient == the_drug, c("study_id", fl), with = FALSE]
  if (!nrow(x)) return(NULL)
  setnames(x, fl, "expo")
  x[, state := fifelse(as.logical(expo), "exposed", "unexposed")]
  e <- ecg_narrow[x[, .(study_id, state)], on = "study_id", nomatch = NULL]
  e[!is.na(state)]
}

# motore "events": distanza dall'ultima somministrazione
classify_events <- function(the_drug, mode = "window",
                            window_h = CFG$exposure_window_h,
                            washout_h = NULL, level = CFG$exposure_level) {
  ad <- admin[ingredient == the_drug, .(subject_id, stay_id, charttime = admin_time)]
  if (!nrow(ad)) return(NULL)
  jk <- if (level == "subject") c("subject_id","charttime")
        else                    c("subject_id","stay_id","charttime")

  e <- copy(ecg_narrow)
  setnames(e, "ecg_time", "charttime")
  e[, row_id := .I]
  setorderv(ad, jk); setorderv(e, jk)
  j <- ad[e, on = jk, roll = TRUE,
          .(row_id = i.row_id, last_adm = x.charttime)]
  e[, last_adm := j[["last_adm"]][match(row_id, j[["row_id"]])]]
  e[, gap_h := as.numeric(difftime(charttime, last_adm, units = "hours"))]
  if (is.null(washout_h)) washout_h <- washout_for(the_drug)

  e[, state := NA_character_]
  e[!is.na(gap_h) & gap_h >= 0 & gap_h <= window_h, state := "exposed"]
  if (mode == "window") {
    e[is.na(state) & (is.na(gap_h) | gap_h > window_h), state := "unexposed"]
  } else if (mode == "washout") {
    e[is.na(state) & (is.na(gap_h) | gap_h > washout_h), state := "unexposed"]
  } else if (mode == "pre_only") {
    e[is.na(state) & is.na(gap_h), state := "unexposed"]
  } else stop("mode non riconosciuto: ", mode)
  setnames(e, "charttime", "ecg_time")
  e[!is.na(state)]
}

classify <- function(the_drug, engine, mode, window_h) {
  if (engine == "flags") {
    if (mode != "window") return(NULL)   # washout impossibile dai booleani
    classify_flags(the_drug, window_h)
  } else classify_events(the_drug, mode, window_h)
}

## ---------------------------------------------------------------------------
## C. CONTRASTO A LIVELLO DI PAZIENTE
## ---------------------------------------------------------------------------

VALUE_COLS <- function() c("qtcf", "hr", "time_from_icu_h", LAB_COLS)

build_contrast <- function(e_state, contrast_def = "all_ecg",
                           level = CFG$contrast_level) {
  vc <- VALUE_COLS()
  miss <- setdiff(vc, names(e_state))
  if (length(miss)) stop("build_contrast: colonne assenti -> ",
                         paste(miss, collapse = ", "))
  e_state <- copy(e_state)
  if (level == "patient") e_state[, stay_id := "ALL"]

  if (contrast_def == "first_ecg") {
    setorder(e_state, subject_id, stay_id, state, ecg_time)
    e_state <- e_state[, .SD[1L], by = .(subject_id, stay_id, state)]
  }

  cnt <- e_state[, .(n_ecg_state = .N), by = .(subject_id, stay_id, state)]
  agg <- e_state[, lapply(.SD, mean, na.rm = TRUE),
                 by = .(subject_id, stay_id, state), .SDcols = vc]
  agg <- agg[cnt, on = c("subject_id","stay_id","state"), nomatch = NULL]

  keep <- c("subject_id","stay_id", vc, "n_ecg_state")
  ex <- agg[state == "exposed",   ..keep]
  un <- agg[state == "unexposed", ..keep]
  if (!nrow(ex) || !nrow(un)) return(NULL)
  setnames(ex, c(vc,"n_ecg_state"), paste0(c(vc,"n_ecg_state"), "_exposed"))
  setnames(un, c(vc,"n_ecg_state"), paste0(c(vc,"n_ecg_state"), "_unexposed"))
  w <- ex[un, on = c("subject_id","stay_id"), nomatch = NULL]   # inner join
  w <- w[is.finite(qtcf_exposed) & is.finite(qtcf_unexposed)]
  if (!nrow(w)) return(NULL)

  # NIENTE get(): w[[nome]] o e' una colonna o non esiste, e lo dice
  for (v in vc) {
    ce <- paste0(v, "_exposed"); cu <- paste0(v, "_unexposed")
    if (!all(c(ce, cu) %in% names(w)))
      stop("build_contrast: colonne assenti dopo il reshape -> ",
           paste(setdiff(c(ce, cu), names(w)), collapse = ", "))
    set(w, j = paste0("d_", v), value = w[[ce]] - w[[cu]])
  }
  set(w, j = "d_log_ecg_count",
      value = log1p(w[["n_ecg_state_exposed"]]) - log1p(w[["n_ecg_state_unexposed"]]))
  set(w, j = "baseline_qtcf", value = w[["qtcf_unexposed"]])

  if (contrast_def == "earliest_stay") {
    setorder(w, subject_id, stay_id)
    w <- w[, .SD[1L], by = subject_id]
  }
  dcols <- c(paste0("d_", vc), "d_log_ecg_count", "baseline_qtcf")
  ns  <- w[, .(n_stays = .N), by = subject_id]
  pat <- w[, lapply(.SD, mean, na.rm = TRUE), by = subject_id, .SDcols = dcols]
  pat[ns, on = "subject_id", nomatch = NULL]
}

COVARS <- function() c("d_time_from_icu_h", "d_hr",
                       paste0("d_", LAB_COLS), "d_log_ecg_count")

## ---------------------------------------------------------------------------
## D. STIMA
## ---------------------------------------------------------------------------

boot_mean_se <- function(x, B = CFG$B_mean) {
  n <- length(x); if (n < 2L) return(NA_real_)
  sd(replicate(B, mean(x[sample.int(n, n, replace = TRUE)])))
}

fit_drug <- function(pat, drug, label, do_boot = TRUE) {
  cv <- COVARS()
  cc <- pat[complete.cases(pat[, c("d_qtcf", cv), with = FALSE])]
  out <- data.table(drug = drug, spec = label,
                    n_informative = nrow(pat), n_complete = nrow(cc))
  if (nrow(cc) < CFG$min_complete_cases) {
    out[, `:=`(estimate = NA_real_, se_hc3 = NA_real_, se_marginal = NA_real_,
               se_boot = NA_real_, se_ratio = NA_real_, p_hc3 = NA_real_,
               p_marginal = NA_real_, degenerate_covars = "", estimable = FALSE)]
    return(out[])
  }
  sds <- cc[, sapply(.SD, sd, na.rm = TRUE), .SDcols = cv]
  degenerate <- names(sds)[!is.finite(sds) | sds < 1e-10]
  cvu <- setdiff(cv, degenerate)
  cc2 <- copy(cc)
  for (v in cvu) set(cc2, j = v, value = cc2[[v]] - mean(cc2[[v]]))

  f <- as.formula(paste("d_qtcf ~",
        if (length(cvu)) paste(cvu, collapse = " + ") else "1"))
  m <- lm(f, data = as.data.frame(cc2))
  V <- hc3_vcov(m)
  b0 <- unname(coef(m)[1]); s_hc3 <- sqrt(V[1,1])
  s_marg <- cc2[, sd(d_qtcf) / sqrt(.N)]
  s_boot <- if (do_boot) boot_mean_se(cc2[["d_qtcf"]]) else NA_real_

  out[, `:=`(estimate = b0, se_hc3 = s_hc3, se_marginal = s_marg,
             se_boot = s_boot, se_ratio = s_hc3 / s_marg,
             ci_lo_hc3  = b0 - 1.96*s_hc3,  ci_hi_hc3  = b0 + 1.96*s_hc3,
             ci_lo_marg = b0 - 1.96*s_marg, ci_hi_marg = b0 + 1.96*s_marg,
             p_hc3      = 2*pnorm(-abs(b0/s_hc3)),
             p_marginal = 2*pnorm(-abs(b0/s_marg)),
             degenerate_covars = paste(degenerate, collapse = "|"),
             estimable = TRUE)]
  out[]
}

## ---------------------------------------------------------------------------
## ESECUZIONE B-D
## ---------------------------------------------------------------------------

msg("\n== B-D. Griglia di specificazioni ==")

SPECS <- rbindlist(lapply(ENGINES, function(en) {
  modes <- if (en == "flags") "window" else c("window","washout","pre_only")
  wins  <- if (en == "flags") c(6,12,24) else CFG$exposure_window_h
  CJ(engine = en, mode = modes, window_h = wins,
     contrast_def = c("all_ecg","first_ecg","earliest_stay"), sorted = FALSE)
}))
if (CFG$contrast_level == "patient") SPECS <- SPECS[contrast_def != "earliest_stay"]
PRIMARY_SPEC <- sprintf("flags/window/%dh/all_ecg", CFG$exposure_window_h)
if (!"flags" %in% ENGINES)
  PRIMARY_SPEC <- sprintf("events/window/%dh/all_ecg", CFG$exposure_window_h)

contrast_cache <- list(); res_list <- list(); fail_list <- list()

run_one <- function(lab, d, en, md, wh, cd) {
  .tb <- NA_character_
  tryCatch(
    withCallingHandlers({
      es <- classify(d, en, md, wh)
      if (is.null(es) || !nrow(es)) NULL else {
        pat <- build_contrast(es, contrast_def = cd)
        if (is.null(pat) || !nrow(pat)) NULL else {
          r <- fit_drug(pat, d, lab, do_boot = identical(lab, PRIMARY_SPEC))
          r[, role := fifelse(d %in% CONTROLS, "negative_control", "candidate")]
          list(pat = pat, r = r)
        }
      }
    }, error = function(e) {
      .tb <<- paste(vapply(sys.calls(),
        function(cc) paste(deparse(cc), collapse = " "), character(1)),
        collapse = "  >>  ")
    }),
    error = function(e) {
      cl <- if (is.null(conditionCall(e))) NA_character_
            else paste(deparse(conditionCall(e)), collapse = " ")
      msg("    ! %s / %s: %s  [in: %s]", lab, d, conditionMessage(e), cl)
      structure(list(err = conditionMessage(e), call = cl, tb = .tb),
                class = "spec_error")
    })
}

for (i in seq_len(nrow(SPECS))) {
  en <- SPECS$engine[i]; md <- SPECS$mode[i]
  wh <- SPECS$window_h[i]; cd <- SPECS$contrast_def[i]
  lab <- sprintf("%s/%s/%dh/%s", en, md, wh, cd)
  msg("  spec %s", lab)
  for (d in ALL_DRUGS) {
    res <- run_one(lab, d, en, md, wh, cd)
    if (is.null(res)) next
    if (inherits(res, "spec_error")) {
      fail_list[[length(fail_list)+1L]] <-
        data.table(spec = lab, drug = d, message = res[["err"]],
                   call = res[["call"]], traceback = res[["tb"]])
      next
    }
    contrast_cache[[paste(lab, d, sep = "|")]] <- res[["pat"]]
    res_list[[length(res_list)+1L]] <- res[["r"]]
  }
}

if (length(fail_list)) {
  FAILS <- rbindlist(fail_list); wr(FAILS, "00_failures")
  msg("\n  %d combinazioni fallite. Espressione del primo errore: %s",
      nrow(FAILS), FAILS[["call"]][1])
  print(FAILS[, .(n = .N), by = message])
}

RES <- rbindlist(res_list, fill = TRUE)
if (!nrow(RES) || !"estimable" %in% names(RES))
  stop("Nessuna stima prodotta: guardare i messaggi qui sopra.")
wr(RES, "01_estimates_by_specification")
if (!RES[estimable == TRUE, .N])
  stop("Tutti i farmaci sotto ", CFG$min_complete_cases, " casi completi: ",
       "controllare la copertura dei laboratori (max_lab_gap_min).")

frag <- RES[estimable == TRUE, .(drug, role, spec, estimate)]
frag <- Reduce(function(a,b) merge(a, b, by = c("drug","role"), all = TRUE),
  lapply(sort(unique(frag$spec)), function(sp) {
    x <- frag[spec == sp, .(drug, role, estimate)]
    setnames(x, "estimate", sp)[]
  }))
wr(frag, "02_fragility_across_specifications")
msg("\n  Fragilita' tra specificazioni:"); print(frag)

msg("\n  SE a confronto (%s):", PRIMARY_SPEC)
print(RES[spec == PRIMARY_SPEC & estimable == TRUE,
          .(drug, role, estimate = round(estimate,2), se_hc3 = round(se_hc3,3),
            se_marginal = round(se_marginal,3), se_boot = round(se_boot,3),
            ratio = round(se_ratio,3))])

## ---------------------------------------------------------------------------
## E. PSEUDO-ESPOSIZIONE TEMPORALE
## ---------------------------------------------------------------------------
# Sugli stessi pazienti, senza alcun farmaco: contrasto "ECG tardivi meno ECG
# precoci" dentro il ricovero. Se riproduce il displacement dei controlli
# negativi, i -2.05 ms sono campionamento temporale, non confondimento.

msg("\n== E. Pseudo-esposizione temporale ==")

pseudo_time_contrast <- function(subjects) {
  e <- ecg_narrow[subject_id %in% subjects]
  if (!nrow(e)) return(NULL)
  grp <- if (CFG$contrast_level == "patient") "subject_id" else c("subject_id","stay_id")
  e <- copy(e)
  e[, med_t := median(ecg_time), by = grp]
  e[, state := fifelse(ecg_time > med_t, "exposed", "unexposed")]
  e <- e[, if (uniqueN(state) == 2L) .SD else NULL, by = grp]
  if (!nrow(e)) return(NULL)
  build_contrast(e, "all_ecg")
}

pseudo_res <- rbindlist(lapply(ALL_DRUGS, function(d) {
  key <- paste(PRIMARY_SPEC, d, sep = "|")
  if (is.null(contrast_cache[[key]])) return(NULL)
  pt <- pseudo_time_contrast(contrast_cache[[key]]$subject_id)
  if (is.null(pt) || !nrow(pt)) return(NULL)
  data.table(drug = d,
             role = fifelse(d %in% CONTROLS, "negative_control", "candidate"),
             n = nrow(pt),
             pseudo_estimate = pt[, mean(d_qtcf, na.rm = TRUE)],
             pseudo_se = pt[, sd(d_qtcf, na.rm = TRUE)/sqrt(.N)])
}), fill = TRUE)

if (nrow(pseudo_res)) {
  pseudo_res <- merge(pseudo_res,
                      RES[spec == PRIMARY_SPEC, .(drug, real_estimate = estimate)],
                      by = "drug", all.x = TRUE)
  wr(pseudo_res, "03_pseudo_time_exposure")
  print(pseudo_res)
  msg("  displacement medio controlli reali: %+.2f ms",
      RES[spec == PRIMARY_SPEC & role == "negative_control",
          mean(estimate, na.rm = TRUE)])
  msg("  displacement medio pseudo-temporale: %+.2f ms",
      pseudo_res[role == "negative_control", mean(pseudo_estimate, na.rm = TRUE)])
}

## ---------------------------------------------------------------------------
## F. STRATIFICAZIONE PER QTcF BASALE
## ---------------------------------------------------------------------------

msg("\n== F. QTcF basale: candidati E controlli ==")

strat_list <- list(); slope_list <- list()
for (d in ALL_DRUGS) {
  pat <- contrast_cache[[paste(PRIMARY_SPEC, d, sep = "|")]]
  if (is.null(pat)) next
  pat <- copy(pat)[is.finite(baseline_qtcf) & is.finite(d_qtcf)]
  if (nrow(pat) < CFG$min_complete_cases) next
  pat[, stratum := fifelse(baseline_qtcf < CFG$baseline_qtcf_cut,
                           sprintf("<%d", CFG$baseline_qtcf_cut),
                           sprintf(">=%d", CFG$baseline_qtcf_cut))]
  s <- pat[, .(n = .N, mean_d = mean(d_qtcf), se = sd(d_qtcf)/sqrt(.N)), by = stratum]
  s[, `:=`(drug = d, role = fifelse(d %in% CONTROLS, "negative_control", "candidate"),
           ci_lo = mean_d - 1.96*se, ci_hi = mean_d + 1.96*se)]
  strat_list[[d]] <- s

  pat[, bl_c := baseline_qtcf - mean(baseline_qtcf)]
  ml <- lm(d_qtcf ~ bl_c, data = as.data.frame(pat))
  Vm <- hc3_vcov(ml)
  sl <- unname(coef(ml)[2]); sl_se <- sqrt(Vm[2,2])
  slope_list[[d]] <- data.table(
    drug = d, role = fifelse(d %in% CONTROLS, "negative_control", "candidate"),
    n = nrow(pat), slope = sl, slope_se = sl_se,
    slope_p = 2*pnorm(-abs(sl/sl_se)),
    intercept_at_mean_baseline = unname(coef(ml)[1]))
}
wr(rbindlist(strat_list, fill = TRUE),  "04_baseline_qtcf_strata")
SLOPES <- rbindlist(slope_list, fill = TRUE); wr(SLOPES, "05_baseline_qtcf_slopes")
if (nrow(SLOPES)) {
  msg("  pendenza di delta su QTcF basale:")
  print(SLOPES[, .(slope_mean = round(mean(slope),3),
                   slope_sd = round(sd(slope),3), n_drugs = .N), by = role])
}

## ---------------------------------------------------------------------------
## G. CALIBRAZIONE
## ---------------------------------------------------------------------------

msg("\n== G. Calibrazione empirica ==")

fit_syserr <- function(y, se) {
  nll <- function(p) {
    v <- exp(p[2])^2 + se^2
    0.5*sum(log(2*pi*v) + (y - p[1])^2/v)
  }
  o <- optim(c(mean(y), log(max(sd(y), 1e-3))), nll, method = "BFGS")
  sg <- exp(o$par[2])
  if (sg < 0.05) warning("fit_syserr: sigma al bordo (", signif(sg,3), ")")
  list(mu = o$par[1], sigma = sg)
}

primary <- RES[spec == PRIMARY_SPEC & estimable == TRUE]
ctrl <- primary[role == "negative_control"]; cand <- primary[role == "candidate"]

calibrate_set <- function(ctrl_dt, cand_dt, set_label, se_col = "se_marginal") {
  if (nrow(ctrl_dt) < 3L || !nrow(cand_dt)) return(NULL)
  f <- fit_syserr(ctrl_dt[["estimate"]], ctrl_dt[[se_col]])
  cd <- copy(cand_dt)
  cd[, `:=`(set = set_label, mu = f[["mu"]], sigma = f[["sigma"]],
            cal_estimate = estimate - f[["mu"]])]
  set(cd, j = "cal_se", value = sqrt(f[["sigma"]]^2 + cd[[se_col]]^2))
  cd[, `:=`(cal_lo = cal_estimate - 1.96*cal_se,
            cal_hi = cal_estimate + 1.96*cal_se,
            cal_p  = 2*pnorm(-abs(cal_estimate/cal_se)))]
  cd[, cal_q := p.adjust(cal_p, method = "BH")]
  cd[, .(set, drug, n_controls = nrow(ctrl_dt), mu, sigma, estimate,
         cal_estimate, cal_lo, cal_hi, cal_p, cal_q)]
}

DROP_CONDITIONAL <- c("famotidine", "pantoprazole")  # CredibleMeds Conditional
calib <- list(calibrate_set(ctrl, cand, "all_controls"),
              calibrate_set(ctrl[!drug %in% DROP_CONDITIONAL], cand,
                            "drop_famotidine_pantoprazole"))
for (d in ctrl$drug)
  calib[[length(calib)+1L]] <- calibrate_set(ctrl[drug != d], cand,
                                             paste0("loo_omit_", d))
CALIB <- rbindlist(calib, fill = TRUE); wr(CALIB, "06_calibration_sets")
if (nrow(CALIB)) {
  msg("  mu e sigma per insieme:")
  print(unique(CALIB[, .(set, n_controls, mu = round(mu,3), sigma = round(sigma,3))]))
}

# bootstrap congiunto: ricampiona i pazienti UNA volta per replica e ricalcola
# tutte le stime prima di rifittare mu e sigma, cosi' la dipendenza fra farmaci
# e l'incertezza della calibrazione entrano nell'intervallo
if (nrow(primary) >= 4L) {
  msg("  bootstrap congiunto (%d repliche)...", CFG$B_calib)
  keys <- setNames(paste(PRIMARY_SPEC, primary$drug, sep = "|"), primary$drug)
  all_subjects <- unique(unlist(lapply(keys, function(k) contrast_cache[[k]]$subject_id)))
  boot_rows <- vector("list", CFG$B_calib)
  for (b in seq_len(CFG$B_calib)) {
    ids <- data.table(subject_id = sample(all_subjects, length(all_subjects), TRUE))
    stat <- lapply(primary$drug, function(d) {
      x <- contrast_cache[[keys[[d]]]][ids, on = "subject_id", nomatch = NULL,
                                       allow.cartesian = TRUE]$d_qtcf
      x <- x[is.finite(x)]
      if (length(x) < 10L) c(NA_real_, NA_real_) else c(mean(x), sd(x)/sqrt(length(x)))
    })
    est <- vapply(stat, `[`, numeric(1), 1L); sev <- vapply(stat, `[`, numeric(1), 2L)
    isc <- primary$drug %in% ctrl$drug & is.finite(est) & is.finite(sev)
    if (sum(isc) < 3L) next
    f <- tryCatch(fit_syserr(est[isc], sev[isc]), error = function(e) NULL)
    if (is.null(f)) next
    boot_rows[[b]] <- data.table(rep = b, drug = primary$drug, est = est,
                                 mu = f[["mu"]], sigma = f[["sigma"]],
                                 cal = est - f[["mu"]])
  }
  BOOT <- rbindlist(boot_rows, fill = TRUE)
  if (nrow(BOOT)) {
    wr(BOOT[is.finite(cal), .(boot_cal_median = median(cal),
        boot_cal_lo = quantile(cal, .025), boot_cal_hi = quantile(cal, .975)),
        by = drug], "07_calibration_bootstrap_ci")
    msg("  mu bootstrap: %.2f (%.2f a %.2f)",
        BOOT[, median(mu, na.rm = TRUE)],
        BOOT[, quantile(mu, .025, na.rm = TRUE)],
        BOOT[, quantile(mu, .975, na.rm = TRUE)])
  }
}

## ---------------------------------------------------------------------------
## H. FDR
## ---------------------------------------------------------------------------

msg("\n== H. FDR: HC3 condizionale vs marginale ==")

fdr <- copy(cand)
fdr[, `:=`(q_hc3 = p.adjust(p_hc3, "BH"), q_marginal = p.adjust(p_marginal, "BH"))]
fdr[, `:=`(sig_hc3 = q_hc3 < 0.05, sig_marginal = q_marginal < 0.05)]
FDR <- fdr[order(p_marginal),
           .(drug, n_complete, estimate = round(estimate,2),
             se_hc3 = round(se_hc3,3), se_marginal = round(se_marginal,3),
             p_hc3 = signif(p_hc3,3), q_hc3 = signif(q_hc3,3),
             p_marginal = signif(p_marginal,3), q_marginal = signif(q_marginal,3),
             sig_hc3, sig_marginal)]
wr(FDR, "08_fdr_hc3_vs_marginal"); print(FDR)
msg("\n  perdono la significativita' FDR con la SE marginale: %s",
    paste(FDR[sig_hc3 & !sig_marginal, drug], collapse = ", "))

## ---------------------------------------------------------------------------

msg("\n== Riepilogo ==")
deg <- unique(RES[nzchar(degenerate_covars), degenerate_covars])
msg("  covariate degeneri: %s", if (length(deg)) paste(deg, collapse = "; ") else "nessuna")
msg("  output in: %s", normalizePath(CFG$out_dir))
writeLines(capture.output(sessionInfo()), file.path(CFG$out_dir, "sessionInfo.txt"))

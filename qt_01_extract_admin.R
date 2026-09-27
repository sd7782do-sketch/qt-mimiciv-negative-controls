###############################################################################
# qt_01_extract_admin.R
#
# Estrae gli EVENTI DI SOMMINISTRAZIONE CON TIMESTAMP da MIMIC-IV.
#
# Serve perche' phase2_ecg_drug_long.parquet contiene solo i booleani
# exposed_6h / exposed_12h / exposed_24h: da quelli non si puo' ricostruire
# quante ore sono passate dall'ultima somministrazione, quindi niente washout
# e niente contrasto pre/post pulito. Qui si recuperano i tempi.
#
# Sorgenti (tutte gia' presenti in C:/Users/sd778/Desktop/mimiciv):
#   emar.csv.gz         somministrazioni registrate (medication, charttime)
#   inputevents.csv.gz  infusioni ICU (itemid, starttime, stay_id)
#   d_items.csv.gz      dizionario itemid -> label
#   icustays.csv.gz     per mappare eMAR (che non ha stay_id) al ricovero ICU
#   qt_feasibility_audit/feas_00_drug_map.csv   ingredient, role, pattern
#
# Output in out_dir:
#   qt_admin_events.parquet     subject_id, hadm_id, stay_id, ingredient,
#                               role, admin_time, source
#   qt_admin_coverage.csv       eventi/pazienti/ricoveri per farmaco
#   qt_emar_event_audit.csv     frequenze di event_txt (per tarare il filtro)
#
# L'SQL e' stato verificato su DuckDB 1.5.5.
###############################################################################

suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(data.table)
})

CFG <- list(
  base_dir = "C:/Users/sd778/Desktop/mimiciv",
  out_dir  = "C:/Users/sd778/Desktop/mimiciv/qt_reanalysis_input",

  f_emar        = "emar.csv.gz",
  f_inputevents = "inputevents.csv.gz",
  f_ditems      = "d_items.csv.gz",
  f_icustays    = "icustays.csv.gz",
  f_drugmap     = "qt_feasibility_audit/feas_00_drug_map.csv",

  # memoria concessa a DuckDB; alzare se la macchina lo consente
  memory_limit = "8GB",
  threads      = max(1L, parallel::detectCores() - 1L)
)

msg <- function(...) cat(sprintf(...), "\n", sep = "")
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

pth <- function(f) {
  p <- file.path(CFG$base_dir, f)
  if (!file.exists(p)) stop("File non trovato: ", p)
  normalizePath(p, winslash = "/")
}

msg("== Estrazione somministrazioni ==")
msg("  base_dir: %s", CFG$base_dir)

con <- dbConnect(duckdb::duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
dbExecute(con, sprintf("SET memory_limit='%s'", CFG$memory_limit))
dbExecute(con, sprintf("SET threads=%d", CFG$threads))

set_var <- function(name, value)
  dbExecute(con, sprintf("SET VARIABLE %s = '%s'", name, value))

set_var("p_drugmap",     pth(CFG$f_drugmap))
set_var("p_icustays",    pth(CFG$f_icustays))
set_var("p_emar",        pth(CFG$f_emar))
set_var("p_ditems",      pth(CFG$f_ditems))
set_var("p_inputevents", pth(CFG$f_inputevents))

SQL <- "
CREATE OR REPLACE TABLE drug_map AS
SELECT lower(trim(ingredient)) AS ingredient,
       lower(trim(role))       AS role,
       lower(trim(pattern))    AS pattern
FROM read_csv_auto(getvariable('p_drugmap'), header = true);

CREATE OR REPLACE TABLE icustays AS
SELECT subject_id, hadm_id, stay_id,
       CAST(intime AS TIMESTAMP)  AS intime,
       CAST(outtime AS TIMESTAMP) AS outtime
FROM read_csv_auto(getvariable('p_icustays'), header = true);

-- audit dei valori di event_txt: guardarlo PRIMA di fidarsi del filtro sotto
CREATE OR REPLACE TABLE emar_event_audit AS
SELECT event_txt, count(*) AS n
FROM read_csv_auto(getvariable('p_emar'), header = true)
GROUP BY 1 ORDER BY n DESC;

-- eMAR: solo eventi in cui il farmaco e' stato davvero somministrato.
-- eMAR non ha stay_id: l'aggancio al ricovero ICU avviene per finestra
-- temporale su icustays, quindi le somministrazioni di reparto restano fuori.
CREATE OR REPLACE TABLE emar_admin AS
WITH e AS (
  SELECT subject_id, hadm_id, emar_id,
         CAST(charttime AS TIMESTAMP) AS admin_time,
         lower(medication) AS medication,
         event_txt
  FROM read_csv_auto(getvariable('p_emar'), header = true)
  WHERE medication IS NOT NULL
    AND charttime IS NOT NULL
    AND (event_txt ILIKE '%administered%'
         OR event_txt IN ('Started','Restarted','Applied'))
    AND event_txt NOT ILIKE '%not given%'
    AND event_txt NOT ILIKE '%held%'
    AND event_txt NOT ILIKE '%stopped%'
)
SELECT DISTINCT
       e.subject_id, e.hadm_id, i.stay_id,
       m.ingredient, m.role, e.admin_time, 'emar' AS source
FROM e
-- I pattern della mappa sono ESPRESSIONI REGOLARI, non stringhe per LIKE:
-- 'acetaminophen|paracetamol' e '(^|[^s])citalopram' con LIKE non matchano
-- mai. regexp_matches interpreta correttamente sia le alternanze sia le
-- sottostringhe semplici usate per gli altri farmaci.
JOIN drug_map m ON regexp_matches(e.medication, m.pattern)
JOIN icustays i
  ON i.subject_id = e.subject_id
 AND (e.hadm_id IS NULL OR i.hadm_id = e.hadm_id)
 AND e.admin_time >= i.intime
 AND e.admin_time <= i.outtime;

-- inputevents: infusioni (amiodarone, propofol, dexmedetomidina...)
CREATE OR REPLACE TABLE inputevents_admin AS
WITH it AS (
  SELECT itemid, lower(label) AS label
  FROM read_csv_auto(getvariable('p_ditems'), header = true)
),
iv AS (
  SELECT subject_id, hadm_id, stay_id, itemid,
         CAST(starttime AS TIMESTAMP) AS admin_time
  FROM read_csv_auto(getvariable('p_inputevents'), header = true)
  WHERE starttime IS NOT NULL
)
SELECT DISTINCT
       iv.subject_id, iv.hadm_id, iv.stay_id,
       m.ingredient, m.role, iv.admin_time, 'inputevents' AS source
FROM iv
JOIN it ON it.itemid = iv.itemid
JOIN drug_map m ON regexp_matches(it.label, m.pattern);

-- unione e deduplica: stessa sostanza, stesso ricovero, stesso minuto = un
-- evento solo (eMAR e inputevents si sovrappongono per le infusioni)
CREATE OR REPLACE TABLE qt_admin_events AS
WITH u AS (
  SELECT * FROM emar_admin
  UNION ALL
  SELECT * FROM inputevents_admin
),
k AS (SELECT *, date_trunc('minute', admin_time) AS minute_key FROM u)
SELECT subject_id, hadm_id, stay_id, ingredient, role,
       min(admin_time) AS admin_time,
       string_agg(DISTINCT source, '+' ORDER BY source) AS source
FROM k
GROUP BY subject_id, hadm_id, stay_id, ingredient, role, minute_key
ORDER BY subject_id, ingredient, admin_time;

CREATE OR REPLACE TABLE qt_admin_coverage AS
SELECT m.ingredient, m.role,
       count(a.admin_time)          AS n_events,
       count(DISTINCT a.subject_id) AS n_subjects,
       count(DISTINCT a.stay_id)    AS n_stays,
       min(a.admin_time)            AS first_event,
       max(a.admin_time)            AS last_event
FROM drug_map m
LEFT JOIN qt_admin_events a ON a.ingredient = m.ingredient
GROUP BY 1,2
ORDER BY n_events DESC;
"

msg("  esecuzione SQL (puo' richiedere diversi minuti su emar.csv.gz)...")
for (stmt in Filter(function(s) nzchar(trimws(s)), strsplit(SQL, ";\n")[[1]]))
  dbExecute(con, stmt)

audit <- as.data.table(dbGetQuery(con, "SELECT * FROM emar_event_audit"))
cover <- as.data.table(dbGetQuery(con, "SELECT * FROM qt_admin_coverage"))

fwrite(audit, file.path(CFG$out_dir, "qt_emar_event_audit.csv"))
fwrite(cover, file.path(CFG$out_dir, "qt_admin_coverage.csv"))

out_parquet <- normalizePath(file.path(CFG$out_dir, "qt_admin_events.parquet"),
                             winslash = "/", mustWork = FALSE)
dbExecute(con, sprintf(
  "COPY qt_admin_events TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD)", out_parquet))

n <- dbGetQuery(con, "SELECT count(*) AS n FROM qt_admin_events")$n

msg("\n  valori di event_txt in eMAR (controllare che il filtro sia giusto):")
print(head(audit, 20))
msg("\n  copertura per farmaco:")
print(cover)
msg("\n  eventi estratti: %s", format(n, big.mark = "."))
msg("  scritto: %s", out_parquet)

msg("\n  ATTENZIONE: se in qt_emar_event_audit.csv compaiono valori di")
msg("  event_txt che indicano somministrazione ma non sono catturati dal")
msg("  filtro (guardare la colonna n), aggiungerli nella clausola WHERE di")
msg("  emar_admin e rieseguire. Il filtro attuale prende gli eventi che")
msg("  contengono 'administered' piu' Started/Restarted/Applied.")

msg("\n  NOTA: i pattern sono trattati come espressioni regolari, come nella")
msg("  mappa originale. Il pattern dell'acetaminofene ('acetaminophen|")
msg("  paracetamol') cattura anche le associazioni con oppioidi")
msg("  (hydrocodone-acetaminophen, oxycodone-acetaminophen). Questo riproduce")
msg("  la pipeline di fase 2, ma per un controllo negativo e' discutibile:")
msg("  controllare in qt_admin_coverage.csv quanto pesano e considerare un")
msg("  pattern piu' stretto come sensibilita'.")

###############################################################################
# qt_04_fix_missing_drugs.R
#
# Acetaminofene e citalopram hanno prodotto ZERO eventi di somministrazione,
# mentre tutti gli altri 22 farmaci sono stati recuperati (senna 67.330 eventi,
# docusato 104.162). Non e' un problema di aggancio al ricovero: e' il match
# sul nome. Questo script trova la causa e propone la correzione.
#
# Tre cause possibili, in ordine di probabilita':
#   1. pattern vuoto o NA nella mappa -> '%' || NULL || '%' restituisce NULL
#      e il join non produce righe
#   2. pattern che non compare in emar.medication con quella grafia
#   3. il farmaco in MIMIC sta in prescriptions/pharmacy ma non in emar
###############################################################################

suppressPackageStartupMessages({library(DBI); library(duckdb); library(data.table)})

CFG <- list(
  base_dir = "C:/Users/sd778/Desktop/mimiciv",
  f_drugmap = "qt_feasibility_audit/feas_00_drug_map.csv",
  f_emar    = "emar.csv.gz",
  drugs     = c("acetaminophen", "citalopram")
)

msg <- function(...) cat(sprintf(...), "\n", sep = "")
pth <- function(f) normalizePath(file.path(CFG$base_dir, f), winslash = "/")

## --- causa 1: la mappa -------------------------------------------------------
map <- fread(pth(CFG$f_drugmap))
msg("== Righe della mappa per i farmaci mancanti ==")
print(map[tolower(trimws(ingredient)) %in% CFG$drugs])
msg("\n  pattern vuoti o NA nell'intera mappa: %d",
    map[is.na(pattern) | !nzchar(trimws(pattern)), .N])
if (map[is.na(pattern) | !nzchar(trimws(pattern)), .N])
  print(map[is.na(pattern) | !nzchar(trimws(pattern))])

## --- causa 2: cosa c'e' davvero in emar --------------------------------------
con <- dbConnect(duckdb::duckdb()); on.exit(dbDisconnect(con, shutdown = TRUE))
dbExecute(con, "SET memory_limit='8GB'")

msg("\n== Stringhe di emar.medication che contengono i nomi ==")
for (d in CFG$drugs) {
  # match largo: primi 6 caratteri del principio attivo
  stem <- substr(d, 1, 6)
  q <- sprintf("
    SELECT lower(medication) AS medication, event_txt, count(*) AS n
    FROM read_csv_auto('%s', header = true)
    WHERE lower(medication) LIKE '%%%s%%'
    GROUP BY 1,2 ORDER BY n DESC LIMIT 15", pth(CFG$f_emar), stem)
  r <- as.data.table(dbGetQuery(con, q))
  msg("\n  -- %s (ricerca su '%s') --", d, stem)
  if (!nrow(r)) {
    msg("     nessuna riga in emar: il farmaco non e' registrato li'.")
    msg("     Verificare prescriptions.csv.gz / pharmacy.csv.gz: la pipeline")
    msg("     di fase 2 potrebbe averlo recuperato da una sorgente diversa,")
    msg("     nel qual caso va aggiunta a qt_01_extract_admin.R.")
  } else {
    print(r)
    msg("     pattern suggerito per la mappa: '%s'",
        tolower(substr(r[1, medication], 1, 10)))
  }
}

## --- causa 3: il farmaco e' altrove ------------------------------------------
msg("\n== Presenza nelle altre sorgenti ==")
for (f in c("prescriptions.csv.gz", "pharmacy.csv.gz")) {
  p <- file.path(CFG$base_dir, f)
  if (!file.exists(p)) { msg("  assente: %s", f); next }
  col <- if (grepl("prescriptions", f)) "drug" else "medication"
  for (d in CFG$drugs) {
    q <- sprintf("SELECT count(*) AS n FROM read_csv_auto('%s', header = true)
                  WHERE lower(%s) LIKE '%%%s%%'",
                 normalizePath(p, winslash = "/"), col, substr(d, 1, 6))
    n <- tryCatch(dbGetQuery(con, q)$n, error = function(e) NA_integer_)
    msg("  %-22s %-14s righe: %s", f, d, format(n, big.mark = ","))
  }
}

msg("\n== Cosa fare ==")
msg("  Se la causa e' 1, correggere feas_00_drug_map.csv e rilanciare qt_01.")
msg("  Se e' 2, sostituire il pattern con quello suggerito sopra.")
msg("  Se e' 3, i due farmaci restano riportabili solo con il motore 'flags':")
msg("  dichiararlo nei Metodi, come gia' scritto nel manoscritto riscritto.")
msg("  L'acetaminofene e' il controllo negativo con il displacement maggiore")
msg("  (-4.05 ms su 2.670 pazienti): vale la pena recuperarlo.")

###############################################################################
# qt_06_lab_gap.R
#
# Produce la distribuzione dell'intervallo fra il prelievo di laboratorio e
# l'elettrocardiogramma a cui e' stato agganciato. Il manoscritto la cita nei
# Limiti, ma qt_02 stampa i quantili a console senza scriverli su file.
#
# Output in out_dir:
#   E_lab_gap_distribution.csv   una riga per analita, piu' la riga complessiva
#   Figure_S_lab_gap.tiff        distribuzione cumulativa per analita
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(DBI); library(duckdb)
})

CFG <- list(
  base_dir = "C:/Users/sd778/Desktop/mimiciv",
  f_phase3 = "qt_feasibility_phase3/phase3_ecg_lab_nearest.parquet",
  out_dir  = "C:/Users/sd778/Desktop/mimiciv/qt_final",
  lab_cols_use = c("potassium","magnesium","calcium..total","creatinine"),
  thresholds_min = c(60, 120, 240, 360)
)

msg <- function(...) cat(sprintf(...), "\n", sep = "")
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

con <- dbConnect(duckdb::duckdb()); on.exit(dbDisconnect(con, shutdown = TRUE))
labs <- as.data.table(dbGetQuery(con, sprintf(
  "SELECT * FROM read_parquet('%s')",
  normalizePath(file.path(CFG$base_dir, CFG$f_phase3), winslash = "/"))))

if (!"lab_gap_min" %in% names(labs))
  stop("phase3 non contiene lab_gap_min: la distanza non e' ricostruibile a valle.")

labs[, label := make.names(tolower(trimws(label)))]
labs <- labs[label %in% CFG$lab_cols_use & is.finite(lab_gap_min)]
labs[, study_id := as.character(study_id)]
# un valore per ECG e analita, il piu' vicino, come nell'analisi
setorder(labs, study_id, label, lab_gap_min)
labs <- unique(labs, by = c("study_id","label"))

pretty_lab <- c(potassium = "Potassium", magnesium = "Magnesium",
                calcium..total = "Total calcium", creatinine = "Creatinine")

summarise_gap <- function(x, nm) {
  qs <- quantile(x, c(.25,.5,.75,.90,.95,.99), na.rm = TRUE)
  out <- data.table(
    analyte = nm, n = length(x),
    median_min = round(qs[["50%"]]),
    iqr_min = sprintf("%.0f\u2013%.0f", qs[["25%"]], qs[["75%"]]),
    p90_min = round(qs[["90%"]]), p95_min = round(qs[["95%"]]),
    p99_min = round(qs[["99%"]]))
  for (th in CFG$thresholds_min)
    set(out, j = sprintf("within_%dmin_pct", th),
        value = round(100 * mean(x <= th, na.rm = TRUE), 1))
  out[]
}

TAB <- rbindlist(c(
  lapply(CFG$lab_cols_use, function(lb)
    summarise_gap(labs[label == lb, lab_gap_min], pretty_lab[[lb]])),
  list(summarise_gap(labs$lab_gap_min, "All analytes"))))

fwrite(TAB, file.path(CFG$out_dir, "E_lab_gap_distribution.csv"))
print(TAB)
msg("\n  scritto: E_lab_gap_distribution.csv")

## --- figura: distribuzione cumulativa ---------------------------------------
d <- copy(labs)[, analyte := pretty_lab[label]]
p <- ggplot(d, aes(lab_gap_min, colour = analyte)) +
  stat_ecdf(geom = "step", linewidth = .8) +
  geom_vline(xintercept = 240, linetype = "22", colour = "grey35") +
  scale_x_continuous(limits = c(0, 720), breaks = seq(0, 720, 120)) +
  scale_y_continuous(labels = function(v) paste0(100*v, "%")) +
  scale_colour_manual(values = c("Potassium" = "#2166AC", "Magnesium" = "#4DAF4A",
                                 "Total calcium" = "#B2182B", "Creatinine" = "#7F7F7F")) +
  labs(x = "Interval between laboratory measurement and electrocardiogram (minutes)",
       y = "Cumulative proportion of measurements") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93"),
        axis.text = element_text(size = 11, colour = "grey15"),
        legend.position = "bottom", legend.title = element_blank(),
        plot.margin = margin(12, 20, 10, 12))
ggsave(file.path(CFG$out_dir, "Figure_S_lab_gap.tiff"), p,
       width = 7.0, height = 5.0, dpi = 600, compression = "lzw", bg = "white")
msg("  scritto: Figure_S_lab_gap.tiff")

## --- frase pronta per i Limiti ----------------------------------------------
a <- TAB[analyte == "All analytes"]
msg("\n  FRASE PER I LIMITI:")
msg("  \"Laboratory values were matched to the nearest available measurement,")
msg("   at a median of %s minutes (interquartile range %s) from the",
    a$median_min, a$iqr_min)
msg("   electrocardiogram; %.1f%% were within four hours (Supplementary Table S9).\"",
    a$within_240min_pct)

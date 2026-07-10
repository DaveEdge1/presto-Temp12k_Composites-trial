#!/usr/bin/env Rscript
# Emit a REAL-ENSEMBLE proxy_ts.json from a chron-repaired slim cache
# (fts_*.rds). Extends the lipd_to_ts.py schema with:
#   age_ensemble    : n_samples x NCOL matrix (real chronology draws)
#   values_ensemble : n_samples x NCOL matrix (real value ensemble)
# Ensembles are column-subsampled to --ncols (deterministic per record) to
# bound the bundled artifact size. build_fts (run_methods.R) consumes these.
# Usage: Rscript emit_realens_json.R --slim fts_cps.rds --out proxy_ts_real.json [--ncols 100]
suppressPackageStartupMessages(library(jsonlite))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(args == f); if (length(i)) args[i + 1] else d }
SLIM <- getarg("--slim"); OUT <- getarg("--out"); NCOL_KEEP <- as.integer(getarg("--ncols", "100"))
s <- readRDS(SLIM); rec <- s$fTS; lat <- s$lat; lon <- s$lon

subsample_cols <- function(m, n, seed) {
  m <- as.matrix(m); nc <- ncol(m)
  if (nc <= n) return(m)
  set.seed(seed)
  m[, sort(sample.int(nc, n)), drop = FALSE]
}
mat_to_rows <- function(m) {           # matrix -> list of row vectors, NaN->null
  lapply(seq_len(nrow(m)), function(i) {
    r <- as.numeric(m[i, ]); r[!is.finite(r)] <- NA; r
  })
}

con <- file(OUT, "w"); writeLines("[", con); first <- TRUE
for (i in seq_along(rec)) {
  t <- rec[[i]]
  ae <- as.matrix(t$ageEnsemble); ve <- as.matrix(t$paleoData_values)
  if (nrow(ae) != nrow(ve) || nrow(ae) < 4) next
  seed <- sum(utf8ToInt(as.character(t$dataSetName %||% paste0("r", i)))) + i
  aeS <- subsample_cols(ae, NCOL_KEEP, seed)
  veS <- subsample_cols(ve, NCOL_KEEP, seed + 1L)
  ageMed <- apply(ae, 1, median, na.rm = TRUE)
  valMed <- apply(ve, 1, median, na.rm = TRUE)
  u <- suppressWarnings(as.numeric(t$paleoData_uncertainty1sd %||% NA))
  obj <- list(
    id = as.character(t$dataSetName %||% paste0("r", i)),
    dataSetName = as.character(t$dataSetName %||% ""),
    age = ifelse(is.finite(ageMed), ageMed, NA),
    values = ifelse(is.finite(valMed), valMed, NA),
    age_ensemble = mat_to_rows(aeS),
    values_ensemble = mat_to_rows(veS),
    lat = as.numeric(lat[i]), lon = as.numeric(lon[i]),
    units = "degc",
    seasonalityGeneral = as.character(t$seasonalityGeneral %||% ""),
    direction = "positive",   # slim cache already pre-flipped negative-direction
    proxy = as.character(t$paleoData_proxy %||% ""),
    uncertainty1sd = if (is.finite(u)) u else NA
  )
  if (!first) writeLines(",", con); first <- FALSE
  writeLines(toJSON(obj, auto_unbox = TRUE, digits = NA, na = "null"), con)
}
writeLines("]", con); close(con)
cat("[emit_realens_json] wrote", OUT, "(<=", NCOL_KEEP, "ens cols)\n")

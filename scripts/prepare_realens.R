#!/usr/bin/env Rscript
# REAL-ENSEMBLE data prep (container step, alternative to lipd_to_ts.py(pickle)).
# Emits proxy_ts.json from a BUNDLED, chron-repaired, column-subsampled slim
# artifact. The pickle carries only single-vector values; the published methods
# used real per-record age+value ensembles (from the 698 lpd files). This
# bundles those (v1.0.0) so the container reproduces the publication.
#
# Method -> PRE-FILTERED record set (each published driver's exact filter, so no
# in-container unit/tag filtering is needed):
#   dcc, gam, scc : temp12kEnsemble + season + degC (779) -> ensemble_dcc.rds
#   cps, paico    : temp12kEnsemble + season       (821) -> ensemble_cpspaico.rds
# All methods use the real temp12kEnsemble VALUE ensembles (real calibration
# realisations). GAM draws the value ensemble but keeps the paper's own Gaussian
# age model (ignores the bundled age_ensemble). SCC needs a multi-column value
# ensemble to composite (a single vector yields an all-NA composite) and its
# cross-cell MEDIAN gridding (run_methods.R) then reproduces the published curve.
#
# Usage: Rscript prepare_realens.R --method <m> --bundle-dir <dir> --out-json <path>
#   bundle-dir must contain: ensemble_dcc.rds, ensemble_cpspaico.rds (each
#   list(fTS, lat, lon), chron-repaired, column-subsampled).
suppressPackageStartupMessages(library(jsonlite))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(args == f); if (length(i)) args[i + 1] else d }
METHOD <- tolower(getarg("--method", "dcc"))
BUNDLE <- getarg("--bundle-dir", "/app/data/realens")
OUT    <- getarg("--out-json")

# emit real-ensemble proxy_ts.json from the method's pre-filtered slim rds
# (dcc/gam/scc share the degC-filtered set; cps/paico share the no-degC set)
rds_name <- if (METHOD %in% c("dcc", "gam", "scc")) "ensemble_dcc.rds" else "ensemble_cpspaico.rds"
src <- file.path(BUNDLE, rds_name)
if (!file.exists(src)) stop("prepare_realens: missing ", src)
s <- readRDS(src); rec <- s$fTS; lat <- s$lat; lon <- s$lon
mat_to_rows <- function(m) lapply(seq_len(nrow(m)), function(i) { r <- as.numeric(m[i, ]); r[!is.finite(r)] <- NA; r })
con <- file(OUT, "w"); writeLines("[", con); first <- TRUE
for (i in seq_along(rec)) {
  t <- rec[[i]]
  ae <- as.matrix(t$ageEnsemble); ve <- as.matrix(t$paleoData_values)
  if (nrow(ae) != nrow(ve) || nrow(ae) < 4) next
  u <- suppressWarnings(as.numeric(t$paleoData_uncertainty1sd %||% NA))
  obj <- list(
    id = as.character(t$dataSetName %||% paste0("r", i)),
    dataSetName = as.character(t$dataSetName %||% ""),
    age = { a <- apply(ae, 1, median, na.rm = TRUE); ifelse(is.finite(a), a, NA) },
    values = { v <- apply(ve, 1, median, na.rm = TRUE); ifelse(is.finite(v), v, NA) },
    age_ensemble = mat_to_rows(ae),
    values_ensemble = mat_to_rows(ve),
    lat = as.numeric(lat[i]), lon = as.numeric(lon[i]),
    units = "degc",
    seasonalityGeneral = as.character(t$seasonalityGeneral %||% ""),
    direction = "positive",     # slim cache already pre-flipped negative-direction
    proxy = as.character(t$paleoData_proxy %||% ""),
    uncertainty1sd = if (is.finite(u)) u else NA
  )
  if (!first) writeLines(",", con); first <- FALSE
  writeLines(toJSON(obj, auto_unbox = TRUE, digits = NA, na = "null"), con)
}
writeLines("]", con); close(con)
cat(sprintf("[prepare_realens] %s: real-ensemble bundle -> %s\n", METHOD, OUT))

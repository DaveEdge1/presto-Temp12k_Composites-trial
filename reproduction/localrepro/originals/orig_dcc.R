#!/usr/bin/env Rscript
# PATCHED COPY of ScientificDataAnalysis/DCC/DCC.R (Kaufman et al. 2020).
# The METHOD is verbatim; the only deviations are environmental, each marked
# [PATCH]:
#   P1 geoChronR bin() monkey-patch (modern geoChronR renamed binvec= to
#      bin.vec=; publication compositeR uses the legacy name)
#   P2 load TS from the shared ts_all.rds cache instead of re-running
#      readLipd+extractTs on the same 698 files (identical data, ~30 min saved)
#   P3 chron-repair (from reproduction/harness/repro.R): modern lipdR 0.6.0
#      mis-pairs some records' ageEnsemble with the measurement table; re-pair
#      from chron<k>model<k>ensemble<k>.csv. Data-loading correction only.
#   P4 no ~/Desktop progress file; cores/nens/out via args; run inside an
#      output dir so committed references are never clobbered.
# Usage: Rscript orig_dcc.R --cache <ts_all.rds> --lpd <lpddir> --out <dir> \
#          [--nens 500] [--ncores 20] [--seed NA]
suppressWarnings(suppressPackageStartupMessages({
  library(geoChronR); library(compositeR); library(purrr); library(magrittr)
  library(foreach); library(doParallel)
}))
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, d = NULL) { i <- which(args == flag); if (length(i)) args[i + 1] else d }
CACHE  <- getarg("--cache"); LPDDIR <- getarg("--lpd")
OUTDIR <- getarg("--out", "out_orig_dcc"); dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
nens   <- as.integer(getarg("--nens", "500"))
NCORES <- as.integer(getarg("--ncores", "20"))
SEED   <- getarg("--seed", NA)
if (!is.na(SEED) && nzchar(SEED) && SEED != "NA") set.seed(as.integer(SEED))

# [PATCH P1] legacy binvec= keyword for modern geoChronR
local({
  orig <- geoChronR::bin
  if (!"binvec" %in% names(formals(orig))) {
    patched <- function(time, values, bin.vec = NULL, binvec = NULL, bin.fun = mean, ...) {
      if (is.null(bin.vec) && !is.null(binvec)) bin.vec <- binvec
      orig(time = time, values = values, bin.vec = bin.vec, bin.fun = bin.fun, ...)
    }
    utils::assignInNamespace("bin", patched, ns = "geoChronR")
  }
})

# [PATCH P2] load database from shared cache (readLipd+extractTs output)
TS <- readRDS(CACHE)

pullTsVariable <- lipdR::pullTsVariable
sg <- pullTsVariable(TS, variable = "interpretation1_seasonalityGeneral")
ic <- pullTsVariable(TS, "paleoData_inCompilation")
u  <- pullTsVariable(TS, "paleoData_units")

# filter by compilation and seasonality  [verbatim from DCC.R]
tu <- which(tolower(ic) == "temp12kensemble" & (tolower(sg) == "annual" | tolower(sg) == "summeronly" | tolower(sg) == "winteronly") & tolower(u) == "degc")
fTS <- TS[tu]
cat(sprintf("[orig_dcc] filtered records: %d\n", length(fTS)))

# [PATCH P3] chron-repair for modern-lipdR ageEnsemble mis-pairing
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
chron_cache <- new.env(parent = emptyenv())
load_chron_ensemble <- function(dsname) {
  if (exists(dsname, envir = chron_cache, inherits = FALSE)) return(get(dsname, envir = chron_cache))
  lpd_path <- file.path(LPDDIR, paste0(dsname, ".lpd"))
  if (!file.exists(lpd_path)) { assign(dsname, NULL, envir = chron_cache); return(NULL) }
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  csvs <- tryCatch(suppressWarnings({
    unzip(lpd_path, exdir = tmp)
    list.files(tmp, pattern = "chron[0-9]+model[0-9]+ensemble[0-9]+\\.csv$",
               recursive = TRUE, full.names = TRUE)
  }), error = function(e) character(0))
  if (length(csvs) == 0) { assign(dsname, NULL, envir = chron_cache); return(NULL) }
  d <- tryCatch(read.csv(csvs[1], header = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(d) || ncol(d) < 2) { assign(dsname, NULL, envir = chron_cache); return(NULL) }
  ens <- suppressWarnings(as.matrix(d[, -1, drop = FALSE])); mode(ens) <- "numeric"
  assign(dsname, ens, envir = chron_cache); ens
}
n_repaired <- 0
for (k in seq_along(fTS)) {
  v <- fTS[[k]]$paleoData_values; ae <- fTS[[k]]$ageEnsemble
  if (is.null(ae) || NROW(as.matrix(v)) == NROW(as.matrix(ae))) next
  new_ae <- load_chron_ensemble(as.character(fTS[[k]]$dataSetName))
  if (!is.null(new_ae) && NROW(new_ae) == NROW(as.matrix(v))) {
    fTS[[k]]$ageEnsemble <- new_ae; n_repaired <- n_repaired + 1
  }
}
good <- vapply(fTS, function(t) {
  v <- t$paleoData_values; ae <- t[["ageEnsemble"]]
  !is.null(ae) && NROW(as.matrix(v)) == NROW(as.matrix(ae)) && NROW(as.matrix(ae)) >= 4
}, logical(1))
cat(sprintf("[orig_dcc] chron-repaired %d; dropping %d mismatched; %d remain\n",
            n_repaired, sum(!good), sum(good)))
fTS <- fTS[good]

# bin the TS  [verbatim]
binvec  <- seq(-50, to = 12050, by = 100)
binAges <- rowMeans(cbind(binvec[-1], binvec[-length(binvec)]))

latbins <- seq(-90, 90, by = 30)
lat <- lipdR::pullTsVariable(fTS, "geo_latitude")

registerDoParallel(NCORES)   # [PATCH P4] was 4

ensOut <- foreach(i = 1:nens) %dopar% {
  scaled <- c()
  for (lb in 1:(length(latbins) - 1)) {
    fi <- which(lat > latbins[lb] & lat <= latbins[lb + 1])
    tc <- compositeEnsembles(fTS[fi], binvec, spread = TRUE, duration = 3000,
                             searchRange = c(0, 7000), gaussianizeInput = FALSE,
                             ageVar = "ageEnsemble", normalizeVariance = FALSE)
    scaled <- cbind(scaled, tc$composite)
  }
  zonalWeights <- sin(latbins[-1] * pi / 180) - sin(latbins[-length(latbins)] * pi / 180)
  zonalWeights <- zonalWeights / sum(zonalWeights)
  scaledDf <- as.data.frame(scaled)
  names(scaledDf) <- stringr::str_c(latbins[-1], " to ", latbins[-length(latbins)])
  scaledDf$year <- binAges
  scaledDf$GlobalMean <- rowSums(t(t(scaled) * zonalWeights))
  return(scaledDf)
}

allLatMeans <- as.matrix(cbind(binAges, map_dfc(ensOut, extract2, "GlobalMean")))
settings <- paste0(nens, "-", length(latbins) - 1, "bands")
readr::write_csv(x = as.data.frame(allLatMeans),
                 file = file.path(OUTDIR, paste0("globalMean", settings, ".csv")),
                 col_names = FALSE)
cat("[orig_dcc] wrote", file.path(OUTDIR, paste0("globalMean", settings, ".csv")), "\n")

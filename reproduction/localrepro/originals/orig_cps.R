#!/usr/bin/env Rscript
# PATCHED COPY of ScientificDataAnalysis/CPS/cps12k.R (Kaufman et al. 2020).
# METHOD verbatim; environmental patches marked:
#   P1 geoChronR bin() legacy-keyword shim + bare-`bin` alias (scaleComposite
#      calls bare bin())
#   P2 TS from shared ts_all.rds cache (identical data to readLipd+extractTs)
#   P3 chron-repair for modern lipdR ageEnsemble mis-pairing (data loading fix)
#   P4 IO/cores/nens/seed via args; targets read from the clone's CPS dir
#   P5 scaleComposite shim: source-identical to compositeR@1e3e0f2e except
#      `if(is.na(scaleWindow))` -> `if(all(is.na(scaleWindow)))` (errors under
#      R >= 4.2 with the length-2 window cps12k.R always passes)
# Usage: Rscript orig_cps.R --cache <ts_all.rds> --lpd <lpddir> \
#          --targets <cps_dir_with_PAGES2k_csvs> --out <dir> \
#          [--nens 500] [--ncores 20] [--seed NA]
suppressWarnings(suppressPackageStartupMessages({
  library(geoChronR); library(compositeR); library(purrr); library(magrittr)
  library(foreach); library(doParallel)
}))
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, d = NULL) { i <- which(args == flag); if (length(i)) args[i + 1] else d }
CACHE   <- getarg("--cache"); LPDDIR <- getarg("--lpd")
TARGDIR <- getarg("--targets")
OUTDIR  <- getarg("--out", "out_orig_cps"); dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
nens    <- as.integer(getarg("--nens", "500"))
NCORES  <- as.integer(getarg("--ncores", "20"))
SEED    <- getarg("--seed", NA)
if (!is.na(SEED) && nzchar(SEED) && SEED != "NA") set.seed(as.integer(SEED))

# [PATCH P1]
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
bin <- geoChronR::bin

# [PATCH P5] scaleComposite: 1e3e0f2e source with only the is.na guard fixed
scaleComposite <- function(composite,binvec,scaleYears,scaleData,scaleWindow = NA,rescale = TRUE,scaleVariance = TRUE){
  if(NCOL(scaleData) > 1){
    d <- bin(scaleYears,values = scaleData[,sample.int(ncol(scaleData),size = 1)],binvec = binvec)
  }else{
    d <- bin(scaleYears,values = scaleData,binvec = binvec)
  }
  if(all(is.na(scaleWindow))){                      # [was: if(is.na(scaleWindow))]
    scaleWindow <- range(scaleYears)
  }
  good <- which(d$x >= min(scaleWindow) & d$x <= max(scaleWindow))
  dv <- d$y[good]
  m <- mean(dv,na.rm = TRUE)
  if(rescale){ m <- 0 }
  s <- sd(dv,na.rm = TRUE)
  compYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))
  swp <- which(compYears >= min(scaleWindow) & compYears <= max(scaleWindow))
  if(scaleVariance){
    scp <- scale(composite,center = mean(composite[swp],na.rm = TRUE),scale  = sd(composite[swp],na.rm = TRUE))
    scaled <- as.matrix(scp)*s+m
  }else{
    scp <- scale(composite,center = mean(composite[swp],na.rm = TRUE),scale  = FALSE)
    scaled <- as.matrix(scp)+m
  }
  return(scaled)
}

# [PATCH P2]
TS <- readRDS(CACHE)
pullTsVariable <- lipdR::pullTsVariable

# filter timeseries  [verbatim from cps12k.R]
sg <- pullTsVariable(TS, variable = "interpretation1_seasonalityGeneral")
ic <- pullTsVariable(TS, "paleoData_inCompilation")
te <- which(tolower(ic) == "temp12kensemble")
gsg <- which(tolower(sg) == "annual" | tolower(sg) == "summeronly" | tolower(sg) == "winteronly")
tu <- intersect(te, gsg)
fTS <- TS[tu]
cat(sprintf("[orig_cps] filtered records: %d\n", length(fTS)))

# [PATCH P3] chron-repair (identical to orig_dcc.R)
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
cat(sprintf("[orig_cps] chron-repaired %d; dropping %d mismatched; %d remain\n",
            n_repaired, sum(!good), sum(good)))
fTS <- fTS[good]

# bin the TS  [verbatim]
binvec  <- seq(-50, to = 12050, by = 100)
binAges <- rowMeans(cbind(binvec[-1], binvec[-length(binvec)]))

latbins <- seq(-90, 90, by = 30)
lat <- pullTsVariable(fTS, "geo_latitude")

# load in scaling data  [P4: from the clone's CPS dir]
targets <- list.files(TARGDIR, pattern = "PAGES2k", full.names = TRUE)
targetsShort <- list.files(TARGDIR, pattern = "PAGES2k", full.names = FALSE)
sw <- 100
targ <- purrr::map(targets, read.csv)
cat(sprintf("[orig_cps] %d PAGES2k targets\n", length(targets)))

registerDoParallel(NCORES)   # [P4] was 4

ensOut <- foreach(i = 1:nens) %dopar% {
  scaled <- c()
  for (lb in 1:(length(latbins) - 1)) {
    fi <- which(lat > latbins[lb] & lat <= latbins[lb + 1])
    tc <- compositeEnsembles(fTS[fi], binvec, spread = TRUE, duration = 3000,
                             searchRange = c(0, 7000), gaussianizeInput = FALSE,
                             ageVar = "ageEnsemble")
    thisTarget <- which(stringr::str_starts(string = targetsShort,
                        paste0(latbins[lb], "to", latbins[lb + 1], "-scaleWindow", sw, "-PAGES2k.csv")))
    if (length(thisTarget) != 1) stop("target matching problem")
    thisScaled <- scaleComposite(composite = tc$composite, binvec = binvec,
                                 scaleYears = 1950 - targ[[thisTarget]][, 1],
                                 scaleData = targ[[thisTarget]][, -1],
                                 scaleWindow = 1950 - c(0, 2000))
    scaled <- cbind(scaled, thisScaled)
  }
  zonalWeights <- sin(latbins[-1] * pi / 180) - sin(latbins[-length(latbins)] * pi / 180)
  zonalWeights <- zonalWeights / sum(zonalWeights)
  scaledDf <- as.data.frame(scaled)
  names(scaledDf) <- stringr::str_c(latbins[-1], " to ", latbins[-length(latbins)])
  scaledDf$year <- binAges
  scaledDf$GlobalMean <- rowSums(t(t(scaled) * zonalWeights))
  return(scaledDf)
}

allLatMeans <- map_dfc(ensOut, extract2, "GlobalMean")
settings <- paste0(nens, "-", length(latbins) - 1, "bands-", sw, "yr2kwindow")
readr::write_csv(x = allLatMeans,
                 file = file.path(OUTDIR, paste0("globalMean", settings, ".csv")))
cat("[orig_cps] wrote", file.path(OUTDIR, paste0("globalMean", settings, ".csv")), "\n")

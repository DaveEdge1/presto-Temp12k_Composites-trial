#!/usr/bin/env Rscript
# STEP 3: run the SHIPPING template (scripts/run_methods.R, f7268c4 compositeR)
# on the REAL-ensemble slim cache, to test whether the template's
# reimplementation of the published engine reproduces the DCC/CPS references.
#
# The template's run_method reimplements compositeEnsembles via
# sampleEnsembleThenBinTs + standardizeMeanIteratively (a DIFFERENT code path
# than the original drivers' direct compositeEnsembles call). Feeding it the
# same real value + age ensembles isolates that reimplementation's fidelity.
#
# Usage: Rscript run_template_realens.R --method dcc|cps --slim <fts.rds> \
#          --refdata <reference_data> --out <csv> [--nens 500] [--ncores 20]
#          [--seed 42] [--age-var ageEnsemble|age] [--run-methods <path>]
suppressWarnings(suppressPackageStartupMessages({
  library(geoChronR); library(compositeR); library(parallel); library(purrr)
}))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, d = NULL) { i <- which(args == flag); if (length(i)) args[i + 1] else d }
METHOD  <- getarg("--method", "dcc")
SLIM    <- getarg("--slim")
REFDIR  <- getarg("--refdata", "reference_data")
OUT     <- getarg("--out", paste0(METHOD, "_global.csv"))
NENS    <- as.integer(getarg("--nens", "500"))
NCORES  <- as.integer(getarg("--ncores", "20"))
SEED    <- as.integer(getarg("--seed", "42"))
AGEVAR  <- getarg("--age-var", "ageEnsemble")
RUNM    <- getarg("--run-methods", "scripts/run_methods.R")
set.seed(SEED)

# legacy geoChronR bin() keyword shim
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

# source the template helpers (strip the main() trailer)
src <- readLines(RUNM)
drop_from <- which(grepl("^if\\s*\\(sys\\.nframe", src))
if (length(drop_from)) src <- src[seq_len(drop_from[1] - 1L)]
eval(parse(text = src), envir = globalenv())

# load real-ensemble slim cache and adapt to the record shape run_method expects
s <- readRDS(SLIM); rec <- s$fTS; lat <- s$lat; lon <- s$lon
fts <- lapply(seq_along(rec), function(i) {
  t <- rec[[i]]
  list(
    dataSetName = as.character(t$dataSetName),
    paleoData_values = as.matrix(t$paleoData_values),      # REAL value ensemble
    ageEnsemble = as.matrix(t$ageEnsemble),                # REAL age ensemble
    age = apply(as.matrix(t$ageEnsemble), 1, median, na.rm = TRUE),  # median fallback
    paleoData_uncertainty1sd = t$paleoData_uncertainty1sd,
    units = "degc",
    seasonalityGeneral = as.character(t$seasonalityGeneral %||% ""),
    proxyType = as.character(t$paleoData_proxy %||% "")
  )
})
cat(sprintf("[template_realens] %s: %d records, age_var=%s\n", METHOD, length(fts), AGEVAR))

LATBINS <- seq(-90, 90, by = 30)
bandIdx <- findInterval(lat, LATBINS, rightmost.closed = TRUE)
bandIdx[bandIdx < 1 | bandIdx > 6 | !is.finite(lat)] <- NA
gridIdx <- rep(NA_integer_, length(fts))

binvec  <- seq(-50, 12050, by = 100)
binAges <- rowMeans(cbind(binvec[-1], binvec[-length(binvec)]))

load_cps_targets <- function(dir) {
  files <- c("-90to-60", "-60to-30", "-30to0", "0to30", "30to60", "60to90")
  lapply(files, function(b) {
    f <- file.path(dir, paste0(b, "-scaleWindow100-PAGES2k.csv"))
    if (!file.exists(f)) return(NULL)
    df <- read.csv(f, check.names = FALSE)
    list(ages = df[[1]], mat = as.matrix(df[, -1, drop = FALSE]))
  })
}
cps_targets <- if (METHOD == "cps") load_cps_targets(file.path(REFDIR, "cps_targets")) else NULL

cfg <- list(ncores = NCORES, seed = SEED, age_var = AGEVAR,
            ref_start = 1800, ref_end = 1900,
            dcc_duration = 3000, cps_duration = 3000)

t0 <- Sys.time()
res <- run_method(METHOD, fts, bandIdx, gridIdx, binvec, binAges, NENS,
                  cps_targets = cps_targets, cfg = cfg)
cat(sprintf("[template_realens] %s done in %.1fs\n", toupper(METHOD),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

g <- data.frame(binAges = binAges, res$global)
names(g) <- c("binAges", paste0("ens", seq_len(ncol(res$global))))
write.csv(g, OUT, row.names = FALSE)
cat("[template_realens] wrote", OUT, "\n")

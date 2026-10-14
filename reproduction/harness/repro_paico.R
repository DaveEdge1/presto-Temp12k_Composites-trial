#!/usr/bin/env Rscript
# Reproduction driver for PaiCo (Pairwise Comparison). Reads the temp12kEnsemble
# slim cache, runs the template's paico.R against the 6 latitudinal bands +
# Neukom 2k CPS target, writes paico_global.csv for compare.py.
suppressWarnings(suppressPackageStartupMessages({
  library(geoChronR); library(compositeR); library(parallel); library(purrr)
}))
# compositeR was written against an older geoChronR API that used `binvec=...`;
# the version we have (1.1.11) renamed it to `bin.vec=`. Monkey-patch
# geoChronR::bin to accept both keywords so compositeR's internal call works.
local({
  orig <- geoChronR::bin
  patched <- function(time, values, bin.vec = NULL, binvec = NULL, bin.fun = mean, ...) {
    if (is.null(bin.vec) && !is.null(binvec)) bin.vec <- binvec
    orig(time = time, values = values, bin.vec = bin.vec, bin.fun = bin.fun, ...)
  }
  utils::assignInNamespace("bin", patched, ns = "geoChronR")
})
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, d = NULL) { i <- which(args == flag); if (length(i)) args[i + 1] else d }
SLIM    <- getarg("--slim",    "/repro/out/fts_dcc.rds")
OUT     <- getarg("--out",     "/repro/out/paico_global.csv")
REFDIR  <- getarg("--refdata", "/app/reference_data")
NENS    <- as.integer(getarg("--nens",   "10"))
NCORES  <- as.integer(getarg("--ncores", "12"))
SEED    <- as.integer(getarg("--seed",   "42"))
set.seed(SEED)

# Source only the PaiCo definitions. run_methods.R has a main() at end that reads
# /results/proxy_ts.json -- inline the few helpers we need to avoid triggering it.
source(getarg("--paico-src", "/app/scripts/paico.R"))

# 30-deg latitudinal bands; area weight = (sin lat_hi - sin lat_lo)/2
N_BANDS <- 6
BAND_WEIGHTS <- {
  e <- seq(-90, 90, by = 30) * pi / 180
  w <- (sin(e[-1]) - sin(e[-length(e)])) / 2
  w / sum(w)
}
area_weight <- function(bandMat) {
  w <- matrix(BAND_WEIGHTS, nrow = nrow(bandMat), ncol = N_BANDS, byrow = TRUE)
  w[!is.finite(bandMat)] <- NA
  num <- rowSums(bandMat * w, na.rm = TRUE)
  den <- rowSums(w, na.rm = TRUE)
  out <- num / den; out[den == 0] <- NA; out
}
apply_reference <- function(ens, binAges, ref_start_ce = 1800, ref_end_ce = 1900,
                            member_ref_bp = NULL) {   # member_ref_bp ignored: full-record default
  ens <- sweep(ens, 2, colMeans(ens, na.rm = TRUE), "-")
  ref_bp <- c(1950 - ref_end_ce, 1950 - ref_start_ce)
  refrows <- which(binAges >= ref_bp[1] & binAges <= ref_bp[2])
  if (length(refrows) == 0) refrows <- which.min(abs(binAges - 100))
  med <- median(apply(ens[refrows, , drop = FALSE], 2, mean, na.rm = TRUE), na.rm = TRUE)
  ens - med
}
load_neukom_targets <- function(dir) {
  if (is.null(dir) || !dir.exists(dir)) return(NULL)
  bands <- c("-90to-60", "-60to-30", "-30to0", "0to30", "30to60", "60to90")
  lapply(bands, function(b) {
    f <- file.path(dir, paste0(b, ".csv"))
    if (!file.exists(f)) return(NULL)
    df <- read.csv(f, check.names = FALSE)
    list(ages = df[["age_bp"]], mat = as.matrix(df[, -1, drop = FALSE]))
  })
}

cat("[repro_paico] loading slim cache ...\n")
s <- readRDS(SLIM)
fts <- s$fTS; lat <- s$lat

# The slim cache stores `age` as a matrix (same shape as ageEnsemble) -- compositeR's
# sampleEnsembleThenBinTs(ageVar="age") considers a non-vector "null" and skips the
# record. Always overwrite `age` with the per-sample median age VECTOR. (PaiCo's
# spread=TRUE then adds ±5%-age noise via the BAM model anyway.)
n_age_fixed <- 0L
for (i in seq_along(fts)) {
  if (!is.null(fts[[i]]$ageEnsemble)) {
    ae <- as.matrix(fts[[i]]$ageEnsemble)
    fts[[i]]$age <- apply(ae, 1, median, na.rm = TRUE)
    n_age_fixed <- n_age_fixed + 1L
  }
}
cat(sprintf("[repro_paico] set age vector (median over ageEnsemble) on %d records\n", n_age_fixed))

LATBINS <- seq(-90, 90, by = 30)
bandIdx <- findInterval(lat, LATBINS, rightmost.closed = TRUE)
bandIdx[bandIdx < 1 | bandIdx > 6 | !is.finite(lat)] <- NA
cat(sprintf("[repro_paico] records: %d, with valid band: %d\n",
            length(fts), sum(!is.na(bandIdx))))

binvec  <- seq(-50, 12050, by = 100)
binAges <- rowMeans(cbind(binvec[-1], binvec[-length(binvec)]))

paico_targets <- load_neukom_targets(file.path(REFDIR, "neukom_targets"))
n_targets_loaded <- sum(!vapply(paico_targets, is.null, logical(1)))
cat(sprintf("[repro_paico] Neukom CPS targets loaded for %d/6 bands\n", n_targets_loaded))

CALIB_LO <- as.numeric(getarg("--calib-lo", "0"))
CALIB_HI <- as.numeric(getarg("--calib-hi", "1000"))   # published last-millennium window (paico.R default)
cfg <- list(ncores = NCORES, paico_reg_param = 100,
            ref_start = 1800, ref_end = 1900,
            seed = SEED,   # per-member RNG stream (matches run_methods.R); also
                           # drives the paper's random per-band target draw
            paico_calib_window = c(CALIB_LO, CALIB_HI))

t0 <- Sys.time()
res <- run_paico(fts, bandIdx, binvec, binAges, NENS,
                 cps_targets = paico_targets,
                 area_weight = area_weight,
                 band_weights = BAND_WEIGHTS,
                 apply_reference = apply_reference,
                 cfg = cfg)
cat(sprintf("[repro_paico] PaiCo done in %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

g <- data.frame(binAges = binAges, res$global)
ne <- ncol(res$global)
names(g) <- c("binAges", paste0("ens", seq_len(ne)))
write.csv(g, OUT, row.names = FALSE)
cat(sprintf("[repro_paico] wrote %s (%d bins x %d members)\n", OUT, nrow(g), ne))

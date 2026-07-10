#!/usr/bin/env Rscript
# Reproduction driver for CPS (Composite Plus Scale). Reads temp12kEnsemble
# slim cache, runs the shippable template's run_method("cps", ...) for the 6
# latitudinal bands, scales each to PAGES2k CPS targets, writes cps_global.csv
# for compare.py. Mirrors repro_paico.R structure.
suppressWarnings(suppressPackageStartupMessages({
  library(geoChronR); library(compositeR); library(parallel); library(purrr)
}))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, d = NULL) { i <- which(args == flag); if (length(i)) args[i + 1] else d }
SLIM    <- getarg("--slim",    "/repro/out/fts_dcc.rds")
OUT     <- getarg("--out",     "/repro/out/cps_global.csv")
REFDIR  <- getarg("--refdata", "/app/reference_data")
NENS    <- as.integer(getarg("--nens",   "10"))
NCORES  <- as.integer(getarg("--ncores", "12"))
SEED    <- as.integer(getarg("--seed",   "42"))
set.seed(SEED)

# Newer geoChronR renamed bin(binvec=) to bin(bin.vec=); compositeR's internal
# call uses the legacy keyword. Monkey-patch so compositeR keeps working.
local({
  orig <- geoChronR::bin
  patched <- function(time, values, bin.vec = NULL, binvec = NULL, bin.fun = mean, ...) {
    if (is.null(bin.vec) && !is.null(binvec)) bin.vec <- binvec
    orig(time = time, values = values, bin.vec = bin.vec, bin.fun = bin.fun, ...)
  }
  utils::assignInNamespace("bin", patched, ns = "geoChronR")
})

# Source the template's run_methods.R helpers we need (don't trigger main()).
# We do this by reading the file and evaluating only the helpers (similar to
# repro_paico.R approach -- inline the few functions to avoid main()).

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
apply_reference <- function(ens, binAges, ref_start_ce = 1800, ref_end_ce = 1900) {
  ens <- sweep(ens, 2, colMeans(ens, na.rm = TRUE), "-")
  ref_bp <- c(1950 - ref_end_ce, 1950 - ref_start_ce)
  refrows <- which(binAges >= ref_bp[1] & binAges <= ref_bp[2])
  if (length(refrows) == 0) refrows <- which.min(abs(binAges - 100))
  med <- median(apply(ens[refrows, , drop = FALSE], 2, mean, na.rm = TRUE), na.rm = TRUE)
  ens - med
}
load_cps_targets <- function(dir) {
  if (is.null(dir) || !dir.exists(dir)) return(NULL)
  files <- c("-90to-60", "-60to-30", "-30to0", "0to30", "30to60", "60to90")
  lapply(files, function(b) {
    f <- file.path(dir, paste0(b, "-scaleWindow100-PAGES2k.csv"))
    if (!file.exists(f)) return(NULL)
    df <- read.csv(f, check.names = FALSE)
    list(ages = df[[1]], mat = as.matrix(df[, -1, drop = FALSE]))
  })
}

# Helpers from run_methods.R: band_composite, scaleCompositeLocal,
# scale_to_target, run_method. run_methods.R has a trailing
# `if (sys.nframe()==0 || identical(environment(), globalenv())) tryCatch(main(),...)`
# that fires under source() and would read /results/proxy_ts.json (missing).
# Strip that block before evaluating.
src <- readLines(getarg("--run-methods", "/app/scripts/run_methods.R"))
# Drop the final block that calls main(). Find the line starting the if-block.
drop_from <- which(grepl("^if\\s*\\(sys\\.nframe", src))
if (length(drop_from)) src <- src[seq_len(drop_from[1] - 1L)]
eval(parse(text = src), envir = globalenv())

cat("[repro_cps] loading slim cache ...\n")
s <- readRDS(SLIM)
fts <- s$fTS; lat <- s$lat

# CPS does NOT require an `age` vector via compositeEnsembles' sampleEnsembleThenBinTs;
# but the slim cache's `age` is a matrix (compositeR rejects). Always overwrite with
# a per-sample median age vector.
for (i in seq_along(fts)) {
  if (!is.null(fts[[i]]$ageEnsemble)) {
    ae <- as.matrix(fts[[i]]$ageEnsemble)
    fts[[i]]$age <- apply(ae, 1, median, na.rm = TRUE)
  }
}
# CPS uses temp12kEnsemble filter (already applied at slim-cache build); does NOT
# require degc-only (the run_method gate is `degc_only <- method %in% c("scc","dcc")`).
# Still: many records have non-degC units; tag .units so the run_method gate works
# even though CPS doesn't gate on units.
for (i in seq_along(fts)) {
  fts[[i]]$units <- "degc"
}

LATBINS <- seq(-90, 90, by = 30)
bandIdx <- findInterval(lat, LATBINS, rightmost.closed = TRUE)
bandIdx[bandIdx < 1 | bandIdx > 6 | !is.finite(lat)] <- NA
cat(sprintf("[repro_cps] records: %d, with valid band: %d\n",
            length(fts), sum(!is.na(bandIdx))))

binvec  <- seq(-50, 12050, by = 100)
binAges <- rowMeans(cbind(binvec[-1], binvec[-length(binvec)]))

cps_targets <- load_cps_targets(file.path(REFDIR, "cps_targets"))
n_targets <- sum(!vapply(cps_targets, is.null, logical(1)))
cat(sprintf("[repro_cps] PAGES2k CPS targets loaded for %d/6 bands\n", n_targets))

# Equal-area grid for SCC reuse — CPS doesn't use it but run_method expects gridIdx
gridIdx <- rep(NA_integer_, length(fts))

cfg <- list(ncores = NCORES, ref_start = 1800, ref_end = 1900,
            cps_duration = 3000, cps_scale_window = c(0, 1000))

t0 <- Sys.time()
res <- run_method("cps", fts, bandIdx, gridIdx, binvec, binAges, NENS,
                  cps_targets = cps_targets, cfg = cfg)
cat(sprintf("[repro_cps] CPS done in %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

g <- data.frame(binAges = binAges, res$global)
ne <- ncol(res$global)
names(g) <- c("binAges", paste0("ens", seq_len(ne)))
write.csv(g, OUT, row.names = FALSE)
cat(sprintf("[repro_cps] wrote %s (%d bins x %d members)\n", OUT, nrow(g), ne))

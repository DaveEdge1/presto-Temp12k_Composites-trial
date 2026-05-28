#!/usr/bin/env Rscript
# Temperature 12k composites -- R methods (DCC, CPS, SCC) + PaiCo dispatch.
#
# Reads proxy_ts.json (from scripts/lipd_to_ts.py), builds a flattened TS list,
# and for each enabled method composites the records within six 30-deg latitude
# bands, then area-weights the bands into a global mean. Each method emits a
# per-method GLOBAL ensemble CSV (binAges + nens columns) and a band-composite
# CSV, consumed downstream by combine_consensus.
#
# Reuses the authors' compositeR engine:
#   - sampleEnsembleThenBinTs : per-member age (BAM) + proxy-unc sampling + binning
#   - standardizeMeanIteratively / standardizeOverInterval : DCC / SCC alignment
#   - scaleComposite : CPS variance scaling to a PAGES2k 2k target
#
# DCC : iterative mean alignment, native degC variance (normalizeVariance=FALSE)
# SCC : align over fixed 3-5 ka window, equal-area gridding, native degC variance
# CPS : z-score (normalizeVariance=TRUE), iterative align, scale to PAGES2k target
# PaiCo: handled in scripts/paico.R (sourced), uses the same binning.

suppressWarnings(suppressPackageStartupMessages({
  library(jsonlite)
  library(yaml)
  library(compositeR)
  library(geoChronR)
  library(purrr)
}))

# ---- six 30-degree bands and their area weights (paper, Methods) ----
LATBINS <- seq(-90, 90, by = 30)                 # 6 bands
BAND_WEIGHTS <- c(0.067, 0.183, 0.25, 0.25, 0.183, 0.067)  # sum = 1.0
N_BANDS <- length(BAND_WEIGHTS)

# ---------------------------------------------------------------------------
# Build the flattened TS list compositeR expects from proxy_ts.json records.
# Pre-applies interpretation direction (flip negative-direction proxies) so we
# don't rely on compositeR's name-based direction detection.
# ---------------------------------------------------------------------------
build_fts <- function(recs) {
  lapply(recs, function(r) {
    vals <- as.numeric(r$values)
    if (identical(tolower(r$direction), "negative")) vals <- vals * -1
    list(
      dataSetName = r$dataSetName,
      age = as.numeric(r$age),
      paleoData_values = vals,
      paleoData_uncertainty1sd = if (is.null(r$uncertainty1sd)) NA_real_ else as.numeric(r$uncertainty1sd),
      lat = as.numeric(r$lat),
      lon = as.numeric(r$lon),
      units = r$units,
      seasonalityGeneral = r$seasonalityGeneral,
      proxyType = r$proxyType
    )
  })
}

band_index <- function(lat) {
  # which 30-deg band each latitude falls in (1..6); NA if out of range
  idx <- findInterval(lat, LATBINS, rightmost.closed = TRUE)
  idx[idx < 1 | idx > N_BANDS | is.na(lat)] <- NA
  idx
}

# nearest equal-area grid cell for each record (great-circle on a sphere)
grid_cell_index <- function(lat, lon, grid) {
  rad <- pi / 180
  glat <- grid$clat * rad; glon <- grid$clon180 * rad
  vapply(seq_along(lat), function(i) {
    if (is.na(lat[i]) || is.na(lon[i])) return(NA_integer_)
    dlon <- glon - lon[i] * rad
    d <- sin(glat) * sin(lat[i] * rad) + cos(glat) * cos(lat[i] * rad) * cos(dlon)
    which.max(pmin(pmax(d, -1), 1))
  }, integer(1))
}

# ---------------------------------------------------------------------------
# Composite the records of one band for one ensemble member.
#   binFunArgs : args forwarded to sampleEnsembleThenBinTs (ar, defaultUnc, ...)
#   stanFun/stanArgs : standardization function + its args
#   gridCells : optional per-record grid-cell ids -> equal-area gridding (SCC)
# Returns the per-bin composite vector (length = #bins).
# ---------------------------------------------------------------------------
band_composite <- function(recs, binvec, binAges, stanFun, stanArgs,
                            binFunArgs = list(), gridCells = NULL) {
  if (length(recs) == 0) return(rep(NA_real_, length(binAges)))

  binMat <- vapply(recs, function(ts) {
    out <- tryCatch(
      do.call(compositeR::sampleEnsembleThenBinTs,
              c(list(ts = ts, binvec = binvec, ageVar = "age",
                     spread = TRUE, alignInterpDirection = FALSE), binFunArgs)),
      error = function(e) rep(NA_real_, length(binAges)))
    if (length(out) != length(binAges)) rep(NA_real_, length(binAges)) else out
  }, numeric(length(binAges)))
  if (is.null(dim(binMat))) binMat <- matrix(binMat, ncol = length(recs))
  binMat[is.nan(binMat)] <- NA

  # need >=2 records with data to standardize
  if (sum(colSums(is.finite(binMat)) > 0) < 2) return(rep(NA_real_, length(binAges)))

  compMat <- tryCatch(
    do.call(stanFun, c(list(ages = binAges, pdm = binMat), stanArgs)),
    error = function(e) binMat)        # fall back to unstandardized on failure
  compMat[!is.finite(compMat)] <- NA

  if (!is.null(gridCells)) {           # equal-area gridding (SCC)
    keep <- which(!is.na(gridCells))
    cm <- compMat[, keep, drop = FALSE]; cg <- gridCells[keep]
    fin <- is.finite(cm); cm0 <- cm; cm0[!fin] <- 0
    # column-group means by grid cell, vectorized via rowsum (records as rows)
    sums <- rowsum(t(cm0), group = cg)
    cnts <- rowsum(t(fin) * 1.0, group = cg)
    cellMat <- t(sums / cnts)            # nbins x ncells; NaN where a cell empty
    cellMat[!is.finite(cellMat)] <- NA
    comp <- rowMeans(cellMat, na.rm = TRUE)
  } else {
    comp <- rowMeans(compMat, na.rm = TRUE)
  }
  comp[is.nan(comp)] <- NA
  comp
}

# area-weight the six band composites into a global mean, renormalizing over
# bands that actually have data at each bin
area_weight <- function(bandMat) {           # bandMat: nbins x 6
  w <- matrix(BAND_WEIGHTS, nrow = nrow(bandMat), ncol = N_BANDS, byrow = TRUE)
  w[!is.finite(bandMat)] <- NA
  num <- rowSums(bandMat * w, na.rm = TRUE)
  den <- rowSums(w, na.rm = TRUE)
  out <- num / den
  out[den == 0] <- NA
  out
}

# subtract each member's full-record mean, then align the ensemble median over
# the reference period (1800-1900 CE -> ~50-150 yr BP) to zero (paper).
apply_reference <- function(ens, binAges, ref_start_ce = 1800, ref_end_ce = 1900) {
  ens <- sweep(ens, 2, colMeans(ens, na.rm = TRUE), "-")
  ref_bp <- c(1950 - ref_end_ce, 1950 - ref_start_ce)   # e.g. c(50,150)
  refrows <- which(binAges >= ref_bp[1] & binAges <= ref_bp[2])
  if (length(refrows) == 0) refrows <- which.min(abs(binAges - 100))
  med <- median(apply(ens[refrows, , drop = FALSE], 2, mean, na.rm = TRUE), na.rm = TRUE)
  ens - med
}

# ---------------------------------------------------------------------------
# Run one method across nens members; returns global ensemble (nbins x nens).
# ---------------------------------------------------------------------------
run_method <- function(method, fts, bandIdx, gridIdx, binvec, binAges, nens,
                       cps_targets = NULL, cfg = list()) {
  degc_only <- method %in% c("scc", "dcc")

  # method-specific standardization + binning settings
  if (method == "dcc") {
    # DCC.R: compositeEnsembles(..., duration=3000, searchRange=c(0,7000),
    #        normalizeVariance=FALSE) with default stan/bin funs.
    stanFun <- compositeR::standardizeMeanIteratively
    stanArgs <- list(duration = cfg$dcc_duration %||% 3000,
                     searchRange = c(0, 7000), normalizeVariance = FALSE)
    binFunArgs <- list(ar = sqrt(0.5))              # AR1 ~0.71 (paper)
  } else if (method == "scc") {
    stanFun <- compositeR::standardizeOverInterval
    stanArgs <- list(interval = c(3000, 5000), normalizeVariance = FALSE)
    binFunArgs <- list(ar = 0)                      # white-noise proxy unc (paper SCC)
  } else if (method == "cps") {
    stanFun <- compositeR::standardizeMeanIteratively
    stanArgs <- list(duration = cfg$cps_duration %||% 3000,
                     searchRange = c(0, 7000), normalizeVariance = TRUE)
    binFunArgs <- list(ar = sqrt(0.5))
  } else {
    stop(paste("run_method: unknown method", method))
  }

  one_member <- function(i) {
    bandMat <- matrix(NA_real_, nrow = length(binAges), ncol = N_BANDS)
    for (b in seq_len(N_BANDS)) {
      sel <- which(bandIdx == b)
      if (degc_only) sel <- sel[vapply(fts[sel], function(x) tolower(x$units) == "degc", logical(1))]
      if (length(sel) < 2) next
      gcells <- if (method == "scc") gridIdx[sel] else NULL
      comp <- band_composite(fts[sel], binvec, binAges, stanFun, stanArgs,
                             binFunArgs = binFunArgs, gridCells = gcells)
      if (method == "cps" && !is.null(cps_targets)) {
        comp <- scale_to_target(comp, binvec, binAges, cps_targets[[b]], cfg)
      }
      bandMat[, b] <- comp
    }
    area_weight(bandMat)
  }

  ncores <- as.integer(cfg$ncores %||% 1)
  if (ncores > 1 && .Platform$OS.type == "unix") {
    cols <- parallel::mclapply(seq_len(nens), one_member, mc.cores = ncores,
                               mc.preschedule = FALSE)
  } else {
    cols <- lapply(seq_len(nens), one_member)
  }
  globalEns <- do.call(cbind, cols)
  apply_reference(globalEns, binAges,
                  cfg$ref_start %||% 1800, cfg$ref_end %||% 1900)
}

# Vendored copy of compositeR's scaleComposite (it is defined in the package but
# NOT exported, so compositeR::scaleComposite errors). Also fixes the original's
# `if (is.na(scaleWindow))` which errors under R >= 4.2 when scaleWindow has
# length 2. Scales the composite's mean+variance over the scaling window to match
# the target's, which restores climate-scale variance to the smoothed composite
# (CPS/PaiCo). Relies on geoChronR::bin being attached.
scaleCompositeLocal <- function(composite, binvec, scaleYears, scaleData,
                                scaleWindow = NA, rescale = TRUE, scaleVariance = TRUE) {
  if (NCOL(scaleData) > 1) {
    d <- bin(scaleYears, values = scaleData[, sample.int(ncol(scaleData), 1)], bin.vec = binvec)
  } else {
    d <- bin(scaleYears, values = scaleData, bin.vec = binvec)
  }
  if (length(scaleWindow) == 1 && is.na(scaleWindow)) scaleWindow <- range(scaleYears)
  good <- which(d$x >= min(scaleWindow) & d$x <= max(scaleWindow))
  dv <- d$y[good]
  m <- if (rescale) 0 else mean(dv, na.rm = TRUE)
  s <- sd(dv, na.rm = TRUE)
  compYears <- rowMeans(cbind(binvec[-1], binvec[-length(binvec)]))
  swp <- which(compYears >= min(scaleWindow) & compYears <= max(scaleWindow))
  cw_sd <- sd(composite[swp], na.rm = TRUE)
  if (!is.finite(cw_sd) || cw_sd == 0) return(composite)
  if (scaleVariance) {
    scp <- scale(composite, center = mean(composite[swp], na.rm = TRUE), scale = cw_sd)
    scaled <- as.matrix(scp) * s + m
  } else {
    scp <- scale(composite, center = mean(composite[swp], na.rm = TRUE), scale = FALSE)
    scaled <- as.matrix(scp) + m
  }
  as.numeric(scaled)
}

# CPS scaling -- replicates cps12k.R exactly:
#   scaleComposite(tc$composite, binvec,
#                  scaleYears  = 1950 - targ[,1],     # target BP -> calendar year
#                  scaleData   = targ[,-1],           # full 2k target ensemble
#                  scaleWindow = 1950 - c(0,2000))    # full 2k window (NOT last-millennium)
# i.e. match the composite's youngest-2k variance to the target's full-2k variance.
# scaleCompositeLocal draws one ensemble column internally when scaleData has >1 col.
scale_to_target <- function(comp, binvec, binAges, target, cfg) {
  if (is.null(target)) return(comp)
  out <- tryCatch(
    scaleCompositeLocal(composite = comp, binvec = binvec,
                        scaleYears = 1950 - target$ages, scaleData = target$mat,
                        scaleWindow = 1950 - c(0, 2000),
                        rescale = TRUE, scaleVariance = TRUE),
    error = function(e) { message("scale_to_target error: ", conditionMessage(e)); comp })
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

load_cps_targets <- function(dir) {
  if (is.null(dir) || !dir.exists(dir)) return(NULL)
  files <- c("-90to-60", "-60to-30", "-30to0", "0to30", "30to60", "60to90")
  lapply(files, function(b) {
    f <- file.path(dir, paste0(b, "-scaleWindow100-PAGES2k.csv"))
    if (!file.exists(f)) return(NULL)
    df <- read.csv(f, check.names = FALSE)
    ages <- df[[1]]
    mat <- as.matrix(df[, -1, drop = FALSE])
    list(ages = ages, mat = mat)
  })
}

# PaiCo scales to the Neukom et al. 2019 CPS 2k reconstruction (the paper's
# targetMedian.CPS), zonal 30-deg band means precomputed from CPS.nc.
# Files: <band>.csv with col "age_bp" (yr BP) + 100 member columns.
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

# ---------------------------------------------------------------------------
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  getarg <- function(flag, default = NULL) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
  }
  ts_path   <- getarg("--ts", "/results/proxy_ts.json")
  cfg_path  <- getarg("--config", "/app/config/user_config.yml")
  refdir    <- getarg("--refdata", "/app/reference_data")
  out_dir   <- getarg("--out-dir", "/results/methods")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cfg <- yaml::read_yaml(cfg_path)
  methods_on <- cfg$methods %||% list(scc = TRUE, dcc = TRUE, cps = TRUE, paico = TRUE, gam = TRUE)
  binc  <- cfg$bin %||% list(start_bp = -50, end_bp = 12050, step = 100)
  nens  <- as.integer(cfg$nens %||% 100)
  refp  <- cfg$reference_period %||% list(start = 1800, end = 1900)

  binvec  <- seq(binc$start_bp, binc$end_bp, by = binc$step)
  binAges <- rowMeans(cbind(binvec[-1], binvec[-length(binvec)]))

  recs <- jsonlite::fromJSON(ts_path, simplifyVector = TRUE,
                             simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
  # keep paper's seasonal subset: annual / summerOnly / winterOnly
  keep_season <- function(r) tolower(r$seasonalityGeneral) %in% c("annual", "summeronly", "winteronly")
  recs <- Filter(keep_season, recs)
  cat(sprintf("[run_methods] %d records after seasonality filter\n", length(recs)))

  fts <- build_fts(recs)
  lat <- vapply(fts, function(x) x$lat, numeric(1))
  lon <- vapply(fts, function(x) x$lon, numeric(1))
  bandIdx <- band_index(lat)

  grid <- read.csv(file.path(refdir, "equal_area_grid_centers.csv"))
  gridIdx <- grid_cell_index(lat, lon, grid)

  cps_targets <- load_cps_targets(file.path(refdir, "cps_targets"))
  # PaiCo scales to the Neukom 2k CPS reconstruction; fall back to PAGES2k if absent.
  paico_targets <- load_neukom_targets(file.path(refdir, "neukom_targets"))
  if (is.null(paico_targets) || all(vapply(paico_targets, is.null, logical(1))))
    paico_targets <- cps_targets

  ncores <- as.integer(cfg$ncores %||% max(1L, parallel::detectCores() - 1L))
  cat(sprintf("[run_methods] using %d core(s)\n", ncores))

  method_cfg <- list(
    dcc_duration = cfg$advanced$dcc_duration,
    cps_duration = cfg$advanced$cps_duration,
    cps_scale_window = cfg$advanced$cps_scale_window,
    paico_reg_param = cfg$advanced$paico_reg_param,
    ncores = ncores,
    ref_start = refp$start, ref_end = refp$end)

  for (m in c("scc", "dcc", "cps")) {
    if (!isTRUE(methods_on[[m]])) next
    cat(sprintf("[run_methods] running %s (nens=%d) ...\n", toupper(m), nens))
    t0 <- Sys.time()
    ens <- run_method(m, fts, bandIdx, gridIdx, binvec, binAges, nens,
                      cps_targets = cps_targets, cfg = method_cfg)
    out <- data.frame(binAges = binAges, ens)
    names(out) <- c("binAges", paste0("ens", seq_len(nens)))
    write.csv(out, file.path(out_dir, paste0(m, "_global.csv")), row.names = FALSE)
    cat(sprintf("[run_methods] %s done in %.1fs -> %s_global.csv\n",
                toupper(m), as.numeric(difftime(Sys.time(), t0, units = "secs")), m))
  }

  # PaiCo (sourced) -- same binning, pairwise-comparison reconstruction
  if (isTRUE(methods_on$paico)) {
    paico_src <- file.path(dirname(getarg("--self", "/app/scripts/run_methods.R")), "paico.R")
    if (!file.exists(paico_src)) paico_src <- "/app/scripts/paico.R"
    if (file.exists(paico_src)) {
      source(paico_src)
      cat("[run_methods] running PaiCo ...\n")
      ens <- run_paico(fts, bandIdx, binvec, binAges, nens,
                       cps_targets = paico_targets, area_weight = area_weight,
                       band_weights = BAND_WEIGHTS, apply_reference = apply_reference,
                       cfg = method_cfg)
      out <- data.frame(binAges = binAges, ens)
      names(out) <- c("binAges", paste0("ens", seq_len(nens)))
      write.csv(out, file.path(out_dir, "paico_global.csv"), row.names = FALSE)
      cat("[run_methods] PaiCo done -> paico_global.csv\n")
    } else {
      cat("[run_methods] paico.R not found; skipping PaiCo\n")
    }
  }

  cat("[run_methods] complete.\n")
}

if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  tryCatch(main(), error = function(e) { message("ERROR: ", conditionMessage(e)); quit(status = 1) })
}

#!/usr/bin/env Rscript
# FAITHFUL-REPRODUCTION HARNESS (not the shippable template).
#
# Reproduces the published Kaufman 2020 per-method GMST by running the authors'
# own compositeR engine on the canonical LiPD files WITH real chronology
# ensembles (ScientificDataAnalysis/lipdFilesWithEnsembles), applying the exact
# published record filter (paleoData_inCompilation == "temp12kEnsemble" &
# seasonalityGeneral in {annual,summerOnly,winterOnly} & degC), with
# ageVar="ageEnsemble" (NOT BAM). Compare output to reference_data/published.
#
# This is reproduction-only. The temp12kEnsemble filter / ensemble-file source
# must NOT leak into the user-facing template (see feedback_reproduction_vs_template).
#
# Usage (in the presto-temp12k container, cwd = / so renv activates):
#   Rscript /repro/repro.R --method dcc --nens 50 --out /repro/out
suppressWarnings(suppressPackageStartupMessages({
  library(lipdR); library(geoChronR); library(compositeR); library(purrr)
}))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
METHOD <- getarg("--method", "dcc")
NENS   <- as.integer(getarg("--nens", "50"))
OUT    <- getarg("--out", "/repro/out")
LPDDIR <- getarg("--lpd", "/repro/T12k/ScientificDataAnalysis/lipdFilesWithEnsembles")
NCORES <- as.integer(getarg("--ncores", "4"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

LATBINS <- seq(-90, 90, by = 30)
ZONAL_W <- sin(LATBINS[-1] * pi / 180) - sin(LATBINS[-length(LATBINS)] * pi / 180)
ZONAL_W <- ZONAL_W / sum(ZONAL_W)                 # area weights, sum=1
N_BANDS <- length(ZONAL_W)
binvec  <- seq(-50, 12050, by = 100)
binAges <- rowMeans(cbind(binvec[-1], binvec[-length(binvec)]))   # 0..12000

# ---- build (and cache) the filtered, cleaned record set -------------------------
# Slim per-method cache so iterating doesn't reload the ~900MB full-TS rds.
slim <- file.path(OUT, paste0("fts_", METHOD, ".rds"))
if (file.exists(slim)) {
  cat("[repro] loading slim record cache:", slim, "\n")
  s <- readRDS(slim); fTS <- s$fTS; lat <- s$lat; lon <- s$lon
  cat("[repro] records:", length(fTS), "\n")
} else {
  ts_cache <- file.path(OUT, "ts_all.rds")
  if (file.exists(ts_cache)) {
    cat("[repro] loading cached TS:", ts_cache, "\n"); TS <- readRDS(ts_cache)
  } else {
    cat("[repro] reading LiPD files from", LPDDIR, "...\n")
    D <- readLipd(LPDDIR); TS <- extractTs(D); saveRDS(TS, ts_cache)
    cat("[repro] cached", length(TS), "TS columns ->", ts_cache, "\n")
  }
  cat("[repro] total TS columns:", length(TS), "\n")

  season <- tolower(vapply(TS, function(t) as.character(t[["interpretation1_seasonalityGeneral"]] %||% NA)[1], character(1)))
  units  <- tolower(vapply(TS, function(t) as.character(t[["paleoData_units"]] %||% NA)[1], character(1)))
  inEns  <- vapply(TS, function(t) any(tolower(as.character(unlist(t[["paleoData_inCompilation"]]))) == "temp12kensemble"), logical(1))

  degc_methods <- c("dcc", "scc", "cps")           # composite methods need degC
  keep <- inEns & season %in% c("annual", "summeronly", "winteronly") &
          (units == "degc" | !(METHOD %in% degc_methods))
  cat(sprintf("[repro] filter: temp12kEnsemble=%d, +season=%d, +degC=%d -> KEEP %d records\n",
              sum(inEns), sum(inEns & season %in% c("annual","summeronly","winteronly")),
              sum(keep & units == "degc"), sum(keep)))
  fTS <- TS[which(keep)]

  # Newer lipdR/geoChronR returns paleoData_values as a value-ensemble matrix and,
  # for a few records, maps the ageEnsemble onto a different-length axis than the
  # values (e.g. val 81 obs vs ageEnsemble 41). compositeEnsembles draws one column
  # from each and requires NROW(values)==NROW(ageEnsemble); one mismatched record
  # aborts the whole band. Drop those inconsistent records.
  good <- vapply(fTS, function(t) {
    v <- t$paleoData_values; ae <- t[["ageEnsemble"]]
    !is.null(ae) && NROW(as.matrix(v)) == NROW(as.matrix(ae)) && NROW(as.matrix(ae)) >= 4
  }, logical(1))
  cat(sprintf("[repro] dropping %d records NROW(values)!=NROW(ageEnsemble); %d remain\n", sum(!good), sum(good)))
  fTS <- fTS[good]

  # Strip each record to the four fields compositeEnsembles/sampleEnsembleThenBinTs
  # actually touches (dataSetName, paleoData_values, ageEnsemble, paleoData_uncertainty1sd)
  # and pre-flip negative-direction proxies (compositeR's str_detect for "_interpDirection"
  # doesn't match this data's "interpretation1_direction" field, so its alignInterpDirection
  # branch is a silent no-op — must be done explicitly). This shrinks the slim cache ~10x
  # and lets fork workers share memory cleanly.
  lat <- vapply(fTS, function(t) suppressWarnings(as.numeric(t$geo_latitude %||% NA)), numeric(1))
  lon <- vapply(fTS, function(t) suppressWarnings(as.numeric(t$geo_longitude %||% NA)), numeric(1))
  dirs <- vapply(fTS, function(t) tolower(as.character(t$interpretation1_direction %||% "positive"))[1], character(1))
  fTS <- lapply(seq_along(fTS), function(i) {
    t <- fTS[[i]]
    v <- t$paleoData_values; if (!is.matrix(v)) v <- matrix(as.numeric(v), ncol = 1)
    if (identical(dirs[i], "negative")) v <- v * -1
    u <- suppressWarnings(as.numeric(t$paleoData_temperature12kUncertainty %||% NA))
    list(dataSetName = as.character(t$dataSetName),
         paleoData_values = v,
         ageEnsemble = as.matrix(t$ageEnsemble),
         paleoData_uncertainty1sd = if (is.finite(u)) u else NULL)
  })
  saveRDS(list(fTS = fTS, lat = lat, lon = lon), slim)
  cat("[repro] wrote slim cache ->", slim, "  (size:", file.info(slim)$size %/% 1e6, "MB)\n")
}

# ---- SCC equal-area gridding setup (per-cell cell-id, computed once) ------------
cell <- NULL
if (METHOD == "scc") {
  g <- read.csv("/repro/equal_area_grid_centers.csv")
  rad <- pi / 180
  cell <- vapply(seq_along(lat), function(i) {
    if (!is.finite(lat[i]) || !is.finite(lon[i])) return(NA_integer_)
    d <- sin(g$clat * rad) * sin(lat[i] * rad) +
         cos(g$clat * rad) * cos(lat[i] * rad) * cos((g$clon180 - lon[i]) * rad)
    which.max(pmin(pmax(d, -1), 1))
  }, integer(1))
  cat(sprintf("[repro] SCC: %d/%d records assigned to %d unique equal-area cells\n",
              sum(!is.na(cell)), length(cell), length(unique(na.omit(cell)))))
}

# ---- method standardization settings (faithful to the published drivers) --------
stan_args <- switch(METHOD,
  dcc = list(duration = 3000, searchRange = c(0, 7000), normalizeVariance = FALSE),
  cps = list(duration = 3000, searchRange = c(0, 7000), normalizeVariance = TRUE),
  scc = NULL,                                                       # SCC has its own path
  stop("method not yet wired in harness: ", METHOD))

# DCC/CPS one_member: compositeEnsembles via the authors' engine.
one_member_compose <- function(m) {
  bandMat <- matrix(NA_real_, nrow = length(binAges), ncol = N_BANDS)
  for (b in seq_len(N_BANDS)) {
    fi <- which(lat > LATBINS[b] & lat <= LATBINS[b + 1])
    if (length(fi) < 2) next
    tc <- tryCatch(
      do.call(compositeEnsembles, c(list(fTS = fTS[fi], binvec = binvec, spread = TRUE,
              gaussianizeInput = FALSE, ageVar = "ageEnsemble", alignInterpDirection = FALSE), stan_args)),
      error = function(e) { message("band ", b, " err: ", conditionMessage(e)); NULL })
    if (!is.null(tc) && !is.null(tc$composite)) bandMat[, b] <- tc$composite
  }
  bandMat
}

# SCC one_member: bin each record -> equal-area grid -> per-cell anomaly vs 3-5ka -> rowMeans.
# Matches SCC_GMST_122719.m's pipeline (gridding BEFORE standardization), using
# real ageEnsemble draws via the authors' binFun.
one_member_scc <- function(m) {
  bandMat <- matrix(NA_real_, nrow = length(binAges), ncol = N_BANDS)
  for (b in seq_len(N_BANDS)) {
    fi <- which(lat > LATBINS[b] & lat <= LATBINS[b + 1])
    if (length(fi) < 2) next
    bm <- vapply(fTS[fi], function(t) {
      out <- tryCatch(
        sampleEnsembleThenBinTs(ts = t, binvec = binvec, ageVar = "ageEnsemble", spread = TRUE,
                                alignInterpDirection = FALSE),
        error = function(e) rep(NA_real_, length(binAges)))
      if (length(out) != length(binAges)) rep(NA_real_, length(binAges)) else out
    }, numeric(length(binAges)))
    if (is.null(dim(bm))) bm <- matrix(bm, ncol = length(fi))
    bm[!is.finite(bm)] <- NA
    cb <- cell[fi]; keep <- which(!is.na(cb) & colSums(is.finite(bm)) > 0)
    if (length(keep) < 2) next
    bm <- bm[, keep, drop = FALSE]; cb <- cb[keep]
    fin <- is.finite(bm); bm0 <- bm; bm0[!fin] <- 0
    sums <- rowsum(t(bm0), group = cb); cnts <- rowsum(t(fin) * 1.0, group = cb)
    cellMat <- t(sums / cnts); cellMat[!is.finite(cellMat)] <- NA
    cellMat <- tryCatch(
      standardizeOverInterval(ages = binAges, pdm = cellMat,
                              interval = c(3000, 5000), normalizeVariance = FALSE),
      error = function(e) cellMat)
    cellMat[!is.finite(cellMat)] <- NA
    bandMat[, b] <- rowMeans(cellMat, na.rm = TRUE)
  }
  bandMat
}

member_fn <- if (METHOD == "scc") one_member_scc else one_member_compose

cat(sprintf("[repro] running %s with nens=%d on %d cores ...\n", toupper(METHOD), NENS, NCORES))
t0 <- Sys.time()
cols <- if (NCORES > 1 && .Platform$OS.type == "unix")
  parallel::mclapply(seq_len(NENS), member_fn, mc.cores = NCORES, mc.preschedule = TRUE) else
  lapply(seq_len(NENS), member_fn)
cat(sprintf("[repro] composited in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

nb <- length(binAges)
cols <- lapply(cols, function(mm) if (is.matrix(mm) && all(dim(mm) == c(nb, N_BANDS))) mm else matrix(NA_real_, nb, N_BANDS))

# area-weight bands -> global per member (renormalize over bands with data)
area_weight <- function(mm) {
  w <- matrix(ZONAL_W, nrow = nb, ncol = N_BANDS, byrow = TRUE); w[!is.finite(mm)] <- NA
  num <- rowSums(mm * w, na.rm = TRUE); den <- rowSums(w, na.rm = TRUE)
  out <- num / den; out[den == 0] <- NA; out
}
glob <- vapply(cols, area_weight, numeric(nb))      # nb x nens

# archival reference convention: per-member remove full-12k mean, then set the
# ensemble median to 0 at 100 BP (NOAA readme).
glob <- sweep(glob, 2, colMeans(glob, na.rm = TRUE), "-")
r100 <- which.min(abs(binAges - 100))
glob <- glob - median(glob[r100, ], na.rm = TRUE)

df <- data.frame(binAges = binAges, glob)
names(df) <- c("binAges", paste0("ens", seq_len(ncol(glob))))
outfile <- file.path(OUT, paste0(METHOD, "_global.csv"))
write.csv(df, outfile, row.names = FALSE)
cat("[repro] wrote", outfile, "\n")

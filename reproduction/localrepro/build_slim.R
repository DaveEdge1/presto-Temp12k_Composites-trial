#!/usr/bin/env Rscript
# Build a method-specific real-ensemble slim cache from ts_all.rds, matching
# each published driver's record filter, with chron-repair. Mirrors
# reproduction/harness/repro.R's filter+repair but standalone and parameterized.
#   dcc : temp12kEnsemble + season + degC   (779; == fts_dcc.rds)
#   cps : temp12kEnsemble + season          (no degC gate, per cps12k.R)
#   paico: temp12kEnsemble + season          (per PaiCo driver)
# Usage: Rscript build_slim.R --method cps --cache ts_all.rds --lpd <dir> --out fts_cps.rds
suppressWarnings(suppressPackageStartupMessages(library(lipdR)))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(args == f); if (length(i)) args[i + 1] else d }
METHOD <- getarg("--method", "cps"); CACHE <- getarg("--cache"); LPDDIR <- getarg("--lpd")
OUT <- getarg("--out")
TS <- readRDS(CACHE)

season_raw <- vapply(TS, function(t) as.character(t[["interpretation1_seasonalityGeneral"]] %||% NA)[1], character(1))
units <- tolower(vapply(TS, function(t) as.character(t[["paleoData_units"]] %||% NA)[1], character(1)))
match_season <- function(s) { if (is.na(s)) return(FALSE); tolower(s) == "annual" || startsWith(s, "summerO") || startsWith(s, "winterO") }
season_ok <- vapply(season_raw, match_season, logical(1))
in_tag <- vapply(TS, function(t) any(tolower(as.character(unlist(t[["paleoData_inCompilation"]]))) == "temp12kensemble"), logical(1))
degc_methods <- c("dcc", "scc")     # NOTE: cps12k.R does NOT gate on degC
keep <- in_tag & season_ok & (units == "degc" | !(METHOD %in% degc_methods))
cat(sprintf("[build_slim] %s: tag=%d +season=%d -> keep %d\n", METHOD, sum(in_tag), sum(in_tag & season_ok), sum(keep)))
fTS <- TS[which(keep)]

# chron-repair (identical to repro.R)
chron_cache <- new.env(parent = emptyenv())
load_chron_ensemble <- function(dsname) {
  if (exists(dsname, envir = chron_cache, inherits = FALSE)) return(get(dsname, envir = chron_cache))
  lpd_path <- file.path(LPDDIR, paste0(dsname, ".lpd"))
  if (!file.exists(lpd_path)) { assign(dsname, NULL, envir = chron_cache); return(NULL) }
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  csvs <- tryCatch(suppressWarnings({ unzip(lpd_path, exdir = tmp)
    list.files(tmp, pattern = "chron[0-9]+model[0-9]+ensemble[0-9]+\\.csv$", recursive = TRUE, full.names = TRUE)
  }), error = function(e) character(0))
  if (length(csvs) == 0) { assign(dsname, NULL, envir = chron_cache); return(NULL) }
  d <- tryCatch(read.csv(csvs[1], header = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(d) || ncol(d) < 2) { assign(dsname, NULL, envir = chron_cache); return(NULL) }
  ens <- suppressWarnings(as.matrix(d[, -1, drop = FALSE])); mode(ens) <- "numeric"
  assign(dsname, ens, envir = chron_cache); ens
}
nrep <- 0
for (k in seq_along(fTS)) {
  v <- fTS[[k]]$paleoData_values; ae <- fTS[[k]]$ageEnsemble
  if (is.null(ae) || NROW(as.matrix(v)) == NROW(as.matrix(ae))) next
  new_ae <- load_chron_ensemble(as.character(fTS[[k]]$dataSetName))
  if (!is.null(new_ae) && NROW(new_ae) == NROW(as.matrix(v))) { fTS[[k]]$ageEnsemble <- new_ae; nrep <- nrep + 1 }
}
good <- vapply(fTS, function(t) { v <- t$paleoData_values; ae <- t[["ageEnsemble"]]
  !is.null(ae) && NROW(as.matrix(v)) == NROW(as.matrix(ae)) && NROW(as.matrix(ae)) >= 4 }, logical(1))
cat(sprintf("[build_slim] repaired %d; drop %d; %d remain\n", nrep, sum(!good), sum(good)))
fTS <- fTS[good]

lat <- vapply(fTS, function(t) suppressWarnings(as.numeric(t$geo_latitude %||% NA)), numeric(1))
lon <- vapply(fTS, function(t) suppressWarnings(as.numeric(t$geo_longitude %||% NA)), numeric(1))
dirs <- vapply(fTS, function(t) tolower(as.character(t$interpretation1_direction %||% "positive"))[1], character(1))
fTS <- lapply(seq_along(fTS), function(i) {
  t <- fTS[[i]]; v <- t$paleoData_values; if (!is.matrix(v)) v <- matrix(as.numeric(v), ncol = 1)
  if (identical(dirs[i], "negative")) v <- v * -1
  u <- suppressWarnings(as.numeric(t$paleoData_temperature12kUncertainty %||% NA))
  list(dataSetName = as.character(t$dataSetName), paleoData_values = v,
       ageEnsemble = as.matrix(t$ageEnsemble),
       paleoData_uncertainty1sd = if (is.finite(u)) u else NULL,
       paleoData_proxy = as.character(t$paleoData_proxy %||% ""),
       seasonalityGeneral = as.character(t$interpretation1_seasonalityGeneral %||% ""))
})
saveRDS(list(fTS = fTS, lat = lat, lon = lon), OUT)
cat("[build_slim] wrote", OUT, "with", length(fTS), "records\n")

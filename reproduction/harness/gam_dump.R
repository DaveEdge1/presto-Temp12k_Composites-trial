#!/usr/bin/env Rscript
# GAM stage 1 (R): read the temp12kEnsemble slim cache, assign each record to a
# 30-deg band + equal-area grid cell, then for each record draw N pairs of
# (ageEnsemble column, value-ensemble column) and write all the sampled (age,
# temp, cell, band) points into one pooled CSV.  Python loads it, fits one
# pygam.LinearGAM per cell, and samples gam.sample(n_draws=nens).
suppressWarnings(suppressPackageStartupMessages({ library(purrr) }))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, d = NULL) { i <- which(args == flag); if (length(i)) args[i + 1] else d }
N_POOL <- as.integer(getarg("--pool", "200"))    # pool draws per record (matches original's 500 cloud size)
OUT    <- getarg("--out", "/repro/out")
SLIM   <- file.path(OUT, "fts_scc.rds")          # reuse SCC slim (has lat+lon)
GRID   <- "/repro/equal_area_grid_centers.csv"
SEED   <- as.integer(getarg("--seed", "42"))

cat("[gam_dump] loading slim cache ...\n")
s <- readRDS(SLIM)
fTS <- s$fTS; lat <- s$lat; lon <- s$lon
cat("[gam_dump] records:", length(fTS), "\n")

g <- read.csv(GRID); rad <- pi / 180
cell <- vapply(seq_along(lat), function(i) {
  if (!is.finite(lat[i]) || !is.finite(lon[i])) return(NA_integer_)
  d <- sin(g$clat * rad) * sin(lat[i] * rad) +
       cos(g$clat * rad) * cos(lat[i] * rad) * cos((g$clon180 - lon[i]) * rad)
  which.max(pmin(pmax(d, -1), 1))
}, integer(1))
LATBINS <- seq(-90, 90, by = 30)
band <- findInterval(lat, LATBINS, rightmost.closed = TRUE)
band[band < 1 | band > 6 | !is.finite(lat)] <- NA
cat(sprintf("[gam_dump] %d records have cell+band; %d unique cells\n",
            sum(!is.na(cell) & !is.na(band)), length(unique(na.omit(cell)))))

set.seed(SEED)
chunks <- vector("list", length(fTS))
for (i in seq_along(fTS)) {
  if (is.na(cell[i]) || is.na(band[i])) next
  r <- fTS[[i]]
  v  <- r$paleoData_values; ae <- r$ageEnsemble
  if (!is.matrix(v)) v  <- matrix(as.numeric(v),  ncol = 1)
  if (!is.matrix(ae)) ae <- matrix(as.numeric(ae), ncol = 1)
  ncv <- NCOL(v); nca <- NCOL(ae)
  ages_k <- ae[, sample.int(nca, N_POOL, replace = TRUE)]
  vals_k <- v[, sample.int(ncv, N_POOL, replace = TRUE)]
  ok <- is.finite(ages_k) & is.finite(vals_k) & ages_k >= -50 & ages_k <= 12050
  if (!any(ok)) next
  chunks[[i]] <- data.frame(
    age  = as.numeric(ages_k[ok]),
    temp = as.numeric(vals_k[ok]),
    cell = cell[i], band = band[i])
}
df <- do.call(rbind, chunks[!vapply(chunks, is.null, logical(1))])
out_path <- file.path(OUT, "gam_pooled.csv")
data.table::fwrite(df, out_path)
cat(sprintf("[gam_dump] wrote %s  (%.1fM rows, %d cells)\n",
            out_path, nrow(df) / 1e6, length(unique(df$cell))))

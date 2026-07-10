#!/usr/bin/env Rscript
# Emit a proxy_ts.json (lipd_to_ts.py schema, single-vector fields) from the
# v1.0.0 Temp12k-tag records in ts_all.rds, for the template GAM (which reads
# r["values"] single vectors + its own age/sigma perturbation). GAM uses the
# Temp12k tag (single vectors), the same data the production pickle carries.
suppressPackageStartupMessages(library(jsonlite))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- which(args == f); if (length(i)) args[i + 1] else d }
CACHE <- getarg("--cache"); OUT <- getarg("--out"); TAG <- tolower(getarg("--tag", "temp12k"))
TS <- readRDS(CACHE)
ic <- tolower(vapply(TS, function(t) as.character(t[["paleoData_inCompilation"]] %||% NA)[1], character(1)))
keep <- which(ic == TAG)
cat(sprintf("[build_proxyts] %d records with tag '%s'\n", length(keep), TAG))
con <- file(OUT, "w"); writeLines("[", con); first <- TRUE
for (k in keep) {
  t <- TS[[k]]
  v <- t$paleoData_values
  vals <- if (NCOL(v) > 1) apply(as.matrix(v), 1, median, na.rm = TRUE) else as.numeric(v)
  age  <- as.numeric(t$age %||% NA)
  if (length(age) != length(vals)) next
  rec <- list(
    id = as.character(t$paleoData_TSid %||% t$dataSetName %||% ""),
    dataSetName = as.character(t$dataSetName %||% ""),
    age = ifelse(is.finite(age), age, NA),
    values = ifelse(is.finite(vals), vals, NA),
    lat = as.numeric(t$geo_latitude %||% t$geo_meanLat %||% NA),
    lon = as.numeric(t$geo_longitude %||% t$geo_meanLon %||% NA),
    units = as.character(t$paleoData_units %||% ""),
    seasonalityGeneral = as.character(t$interpretation1_seasonalityGeneral %||% ""),
    direction = as.character(t$interpretation1_direction %||% ""),
    proxy = as.character(t$paleoData_proxy %||% t$paleoData_proxyGeneral %||% "")
  )
  if (!first) writeLines(",", con); first <- FALSE
  writeLines(toJSON(rec, auto_unbox = TRUE, digits = NA, na = "null"), con)
}
writeLines("]", con); close(con)
cat("[build_proxyts] wrote", OUT, "\n")

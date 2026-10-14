#!/usr/bin/env Rscript
# Export the tagged TS entries from ts_all.rds to ndjson for scipy->TS.mat.
# Only the fields the published MATLAB drivers touch (SCC_GMST_122719.m,
# PaiCo*). Values/age kept as-is (Temp12k entries are single vectors;
# temp12kEnsemble entries are matrices).
suppressPackageStartupMessages(library(jsonlite))
`%||%` <- function(a,b) if (is.null(a)||length(a)==0||all(is.na(a))) b else a
args <- commandArgs(trailingOnly = TRUE)
CACHE <- args[1]; OUT <- args[2]
TS <- readRDS(CACHE)
ic <- tolower(vapply(TS, function(t) as.character(t[["paleoData_inCompilation"]] %||% NA)[1], character(1)))
keep <- which(ic %in% c("temp12k", "temp12kensemble"))
cat("[export_ts_mat] exporting", length(keep), "tagged TS entries\n")
con <- file(OUT, "w")
for (k in keep) {
  t <- TS[[k]]
  rec <- list(
    dataSetName = as.character(t$dataSetName %||% ""),
    geo_latitude = as.numeric(t$geo_latitude %||% NA),
    geo_meanLat = as.numeric(t$geo_meanLat %||% t$geo_latitude %||% NA),
    geo_meanLon = as.numeric(t$geo_meanLon %||% t$geo_longitude %||% NA),
    paleoData_units = as.character(t$paleoData_units %||% ""),
    paleoData_inCompilation = as.character(t$paleoData_inCompilation %||% "")[1],
    interpretation1_seasonalityGeneral = as.character(t$interpretation1_seasonalityGeneral %||% ""),
    paleoData_temperature12kUncertainty = as.character(t$paleoData_temperature12kUncertainty %||% "NA")[1],
    age = as.numeric(t$age),
    values = if (NCOL(t$paleoData_values) > 1) apply(as.matrix(t$paleoData_values), 1, median, na.rm = TRUE) else as.numeric(t$paleoData_values)
  )
  writeLines(toJSON(rec, digits = NA, auto_unbox = TRUE, na = "null"), con)
}
close(con)
cat("[export_ts_mat] wrote", OUT, "\n")

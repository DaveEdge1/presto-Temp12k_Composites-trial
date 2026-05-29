suppressWarnings(suppressPackageStartupMessages({ library(lipdR); library(geoChronR); library(compositeR); library(purrr) }))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
TS <- readRDS("/repro/out/ts_all.rds")
cat("TS cols:", length(TS), "\n")

inEns <- vapply(TS, function(t) any(tolower(as.character(unlist(t[["paleoData_inCompilation"]]))) == "temp12kensemble"), logical(1))
season <- tolower(vapply(TS, function(t) as.character(t[["interpretation1_seasonalityGeneral"]] %||% NA)[1], character(1)))
units  <- tolower(vapply(TS, function(t) as.character(t[["paleoData_units"]] %||% NA)[1], character(1)))
keep <- inEns & season %in% c("annual","summeronly","winteronly") & units == "degc"
fTS <- TS[which(keep)]
cat("kept:", length(fTS), "\n")

lat <- vapply(fTS, function(t) suppressWarnings(as.numeric(t$geo_latitude %||% NA)), numeric(1))
cat("lat non-NA:", sum(is.finite(lat)), "/", length(lat), " range:", toString(round(range(lat, na.rm=TRUE),1)), "\n")

# scan value/age shapes across all kept records
scan <- t(vapply(fTS, function(t) {
  v <- t$paleoData_values; a <- t$age
  c(nval = length(v), nage = length(a),
    is_mat = as.numeric(is.matrix(v)),
    vcols = NCOL(v),
    mult = as.numeric(length(a) > 0 && length(v) %% length(a) == 0 && length(v) != length(a)))
}, numeric(5)))
cat("\n-- value-shape scan over", nrow(scan), "kept records --\n")
cat("len(values)==len(age):", sum(scan[,"nval"]==scan[,"nage"]), "\n")
cat("paleoData_values is.matrix:", sum(scan[,"is_mat"]), "\n")
cat("len(values) multiple of len(age) (not equal):", sum(scan[,"mult"]), "\n")
mm <- which(scan[,"nval"]!=scan[,"nage"] & scan[,"is_mat"]==0)
cat("MISMATCH (flat, len!=age):", length(mm), "\n")
if (length(mm)) {
  i <- mm[1]; r <- fTS[[i]]
  cat("  example mismatch idx", i, "nval", length(r$paleoData_values), "nage", length(r$age),
      "ratio", length(r$paleoData_values)/length(r$age), "\n")
  cat("  is.matrix(values):", is.matrix(r$paleoData_values), " class:", class(r$paleoData_values), "\n")
  cat("  variableName:", toString(r$paleoData_variableName), " units:", toString(r$paleoData_units), "\n")
}

# band 4 (0-30N) records, try compositeEnsembles WITHOUT tryCatch
LATBINS <- seq(-90,90,by=30); binvec <- seq(-50,12050,by=100)
b <- 4
fi <- which(lat > LATBINS[b] & lat <= LATBINS[b+1])
cat("\nband", b, "(", LATBINS[b], "to", LATBINS[b+1], ") n records:", length(fi), "\n")
cat("\n== per-record dims in band", b, "(nrow(val) vs nrow(ageEns) vs len(age)) ==\n")
bad <- 0
for (j in fi) {
  r <- fTS[[j]]
  dv <- dim(as.matrix(r$paleoData_values)); de <- dim(as.matrix(r[["ageEnsemble"]])); la <- length(r$age)
  mism <- dv[1] != de[1]
  if (mism) bad <- bad + 1
  if (mism || j %in% fi[1:3]) cat(sprintf("  idx %d: val %s, ageEns %s, age %d %s\n",
      j, toString(dv), toString(de), la, ifelse(mism, "<-- NROW MISMATCH", "")))
}
cat("records in band with NROW(val)!=NROW(ageEns):", bad, "/", length(fi), "\n")

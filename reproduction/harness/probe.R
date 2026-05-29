suppressWarnings(suppressPackageStartupMessages({
  library(lipdR); library(geoChronR)
}))
D <- readLipd("/repro/probe_subset")
TS <- extractTs(D)
cat("n TS records:", length(TS), "\n")

units <- vapply(TS, function(t) toString(t$paleoData_units), character(1))
cat("\nunits table:\n"); print(table(units))

# find a temperature record (degC)
deg <- which(tolower(units) == "degc")
cat("\ndegC records:", length(deg), "\n")
if (length(deg) == 0) {
  # maybe variableName-based
  vn <- vapply(TS, function(t) toString(t$paleoData_variableName), character(1))
  cat("variableName table (first 20):\n"); print(head(sort(table(vn), decreasing = TRUE), 20))
}

idx <- if (length(deg)) deg[1] else 1
r <- TS[[idx]]
cat("\n== FULL field names of record", idx, "==\n")
print(sort(names(r)))

cat("\n== candidate metadata values on this record ==\n")
for (f in c("paleoData_inCompilation","paleoData_inCompilationBeta",
            "interpretation1_seasonalityGeneral","interpretation1_seasonality",
            "interpretation1_variable","interpretation1_direction",
            "geo_meanLat","geo_latitude","geo_meanLon","geo_longitude",
            "paleoData_variableName","paleoData_units",
            "paleoData_temperature12kUncertainty","archiveType")) {
  cat(sprintf("  %-38s = %s\n", f, toString(r[[f]])))
}

# Try pullTsVariable like the original DCC.R
cat("\n== pullTsVariable probes (original DCC.R fields) ==\n")
for (v in c("paleoData_inCompilation","interpretation1_seasonalityGeneral","paleoData_units","geo_meanLat")) {
  out <- tryCatch(pullTsVariable(TS, v), error = function(e) paste("ERR:", conditionMessage(e)))
  cat(v, "-> non-empty:", sum(nzchar(as.character(out)) & !is.na(out)), "/", length(out),
      " sample:", toString(utils::head(unique(out[nzchar(as.character(out))]), 4)), "\n")
}

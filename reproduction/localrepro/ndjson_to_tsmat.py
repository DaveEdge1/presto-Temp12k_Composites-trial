#!/usr/bin/env python3
"""Convert ts_tagged.ndjson -> TS.mat (1xN MATLAB struct array) for the
published SCC driver. Field names must match SCC_GMST_122719.m exactly."""
import json, sys
import numpy as np
from scipy.io import savemat

src, dst = sys.argv[1], sys.argv[2]
recs = [json.loads(l) for l in open(src)]
n = len(recs)
print(f"[ndjson_to_tsmat] {n} records")

def col(v):
    a = np.asarray(v, dtype=float).reshape(-1, 1)
    return a

fields = {}
fields["dataSetName"] = np.empty((1, n), dtype=object)
fields["geo_latitude"] = np.empty((1, n), dtype=object)
fields["geo_meanLat"] = np.empty((1, n), dtype=object)
fields["geo_meanLon"] = np.empty((1, n), dtype=object)
fields["paleoData_units"] = np.empty((1, n), dtype=object)
fields["paleoData_inCompilation"] = np.empty((1, n), dtype=object)
fields["interpretation1_seasonalityGeneral"] = np.empty((1, n), dtype=object)
fields["paleoData_temperature12kUncertainty"] = np.empty((1, n), dtype=object)
fields["age"] = np.empty((1, n), dtype=object)
fields["paleoData_values"] = np.empty((1, n), dtype=object)

for i, r in enumerate(recs):
    fields["dataSetName"][0, i] = r["dataSetName"]
    fields["geo_latitude"][0, i] = float(r["geo_latitude"]) if r["geo_latitude"] is not None else np.nan
    fields["geo_meanLat"][0, i] = float(r["geo_meanLat"]) if r["geo_meanLat"] is not None else np.nan
    fields["geo_meanLon"][0, i] = float(r["geo_meanLon"]) if r["geo_meanLon"] is not None else np.nan
    fields["paleoData_units"][0, i] = r["paleoData_units"]
    fields["paleoData_inCompilation"][0, i] = r["paleoData_inCompilation"]
    fields["interpretation1_seasonalityGeneral"][0, i] = r["interpretation1_seasonalityGeneral"]
    # MATLAB code does strncmp(er,'NA',2); keep as char when NA, double otherwise
    unc = r["paleoData_temperature12kUncertainty"]
    try:
        fields["paleoData_temperature12kUncertainty"][0, i] = float(unc)
    except (TypeError, ValueError):
        fields["paleoData_temperature12kUncertainty"][0, i] = "NA"
    fields["age"][0, i] = col([np.nan if v is None else v for v in r["age"]])
    fields["paleoData_values"][0, i] = col([np.nan if v is None else v for v in r["values"]])

# Build a MATLAB struct array: savemat maps a numpy structured/object approach;
# simplest reliable route is one cell array per field + cell2struct in MATLAB,
# but savemat CAN write a struct array from a dict-of-object-arrays via
# np.core.records. Instead: save fields separately; the MATLAB wrapper
# assembles TS = struct('age', agec, ...) which builds a 1xN struct array
# from same-sized cell arrays.
savemat(dst, {f"f_{k}": v for k, v in fields.items()}, long_field_names=True, do_compression=True)
print(f"[ndjson_to_tsmat] wrote {dst}")

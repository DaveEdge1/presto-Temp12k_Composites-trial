# GAM port references (for NEXT TASK 1 in ../HANDOFF.md)

- `gam_ensemble_ORIGINAL.py` — the published GAM_frozen engine (nickmckay/
  Temperature12k GAM_frozen/scripts). The algorithm to port: _compute_anomaly
  (modern-window ensemble alignment, L486+) + _predict_gam 0-insertion at -35 BP
  (L236). Reference curve to match: reference_data/published/gam_published.csv.
- `gam_modern_BROKEN_attempt.py` — earlier failed crude modern-anchor (scalar
  worldclim subtraction, aligned unions never anchored -> stayed absolute degC ->
  maxD 10.8). The bug: it worldclim-anchored only solo/no-overlap records, not
  the aligned-union base (left at offset 0.0), mixing absolute unions with
  anomaly solos.
- `gam_modern_FAITHFUL_attempt.py` — the CORRECT faithful port (2026-07-03):
  iterative-union alignment + ONE modern-window (-50..-20 BP) reference per union
  (worldclim fallback) + 0-insertion at -35 BP + keep-all-cells. Env knobs
  (GAM_ANCHOR_MODE, GAM_WC_FALLBACK, GAM_ANCHOR_FRAC) parameterize the variants
  tested. RESULT: the faithful modern-anchor is WORSE than the 3-5 ka baseline
  (see VALIDATION.md "GAM fix attempts" for the table). The 0-insertion adds a
  warm bias; the worldclim-absolute anchoring inflates amplitude. Conclusion:
  anchoring is NOT the dominant GAM residual; the input-data reimplementation
  (synthetic single-vector + sigma noise vs the original's real value ensembles)
  and pygam 0.8->0.12 dominate. **UPDATE: GAM was subsequently FIXED (maxD
  0.259 -> 0.155) by feeding the real per-record VALUE ensembles + a 0-insertion
  at -35 BP into shipping gam_method.py (NOT the modern anchor). See VALIDATION.md
  "GAM FIXED via real VALUE ensembles".** This modern-anchor attempt remains here
  only as the record of the refuted approach.

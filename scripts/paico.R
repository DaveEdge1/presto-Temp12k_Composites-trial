# PaiCo -- Pairwise Comparison reconstruction.
#
# Faithful port of Sami Hanhijarvi's PaiCo solver as archived in
# nickmckay/Temperature12k/ScientificDataAnalysis/PaiCo/Sami/paico/
# (paico.m + private/{inferProxies,selectProxies,optimize,calibrate}.m and the
# C++ core paico_mex.cpp). The C++ is only a performance optimisation for
# sub-annual / overlapping sample windows; once every record is binned onto the
# common target grid (as here), each pairwise comparison is between two grid
# bins, so the comparison matrix A and win/loss counts can be built directly and
# the reconstruction is the maximum-likelihood signal from optimize.m.
#
# Algorithm (per latitude band, per ensemble member):
#   1. bin each record (sampleEnsembleThenBinTs: BAM age + proxy-unc sampling),
#      flip negative-interpDirection records (PaiCoByLat.m).
#   2. for every ordered bin pair (i>j), across records, tally
#        cpos = #(value_i > value_j),  cneg = #(value_i < value_j)
#      and build A with row = e_i - e_j  (z = A f = f_i - f_j).
#   3. Newton-Raphson MLE for the signal f maximising
#        sum( cpos*log Phi(z) + cneg*log(1-Phi(z)) ) - (1/regcov) ||f||^2
#      (probit pairwise-comparison likelihood; optimize.m, adaptive
#       Levenberg-Marquardt damping, heuristic start = row-sum of A).
#   4. calibrate (calibrate.m): mean-variance match f to the 2k target over the
#      overlap window, giving degC. NB the paper's exact target (targetMedian.CPS
#      / Neukom 2k field recon) is not archived; we calibrate to the bundled
#      PAGES2k 2k composite, the closest available 2k target.
# then area-weight the six bands (sin-based zonal weights) and reference.

# ---- pairwise comparison counts + matrix A from a binned record matrix --------
.paico_counts <- function(binMat) {
  nbins <- nrow(binMat)
  Cpos <- matrix(0, nbins, nbins)
  Cneg <- matrix(0, nbins, nbins)
  for (k in seq_len(ncol(binMat))) {
    x <- binMat[, k]
    obs <- which(is.finite(x))
    if (length(obs) < 2) next
    xo <- x[obs]
    dif <- outer(xo, xo, "-")             # dif[a,b] = x_a - x_b
    ut <- which(upper.tri(dif), arr.ind = TRUE)   # a<b in obs-index space
    ia <- obs[ut[, 1]]; ib <- obs[ut[, 2]]; dv <- dif[upper.tri(dif)]
    pos <- dv > 0; neg <- dv < 0
    if (any(pos)) { idx <- cbind(ia[pos], ib[pos]); Cpos[idx] <- Cpos[idx] + 1 }
    if (any(neg)) { idx <- cbind(ia[neg], ib[neg]); Cneg[idx] <- Cneg[idx] + 1 }
  }
  pairs <- which((Cpos + Cneg) > 0, arr.ind = TRUE)   # i<j (upper tri)
  if (nrow(pairs) < 2) return(NULL)
  np <- nrow(pairs)
  A <- matrix(0, np, nbins)
  A[cbind(seq_len(np), pairs[, 1])] <- 1      # +f_i
  A[cbind(seq_len(np), pairs[, 2])] <- -1     # -f_j
  list(A = A, cpos = Cpos[pairs], cneg = Cneg[pairs])
}

# ---- Newton-Raphson MLE (faithful port of optimize.m) -------------------------
.paico_optimize <- function(A, cpos, cneg, regcov_scalar = 1,
                            errorTolerance = 1e-8, maxIters = 1e4L) {
  n <- ncol(A)
  regcov <- diag(n) / regcov_scalar           # paico.m: eye(n)/regcov
  eps <- .Machine$double.eps
  total <- cpos + cneg

  f <- as.numeric(((cpos - cneg) / total) %*% A)   # heuristic start
  s <- stats::sd(f); if (!is.finite(s) || s == 0) s <- 1
  f <- f / s

  z <- as.numeric(A %*% f)
  cdf <- pnorm(z); cdf[cdf == 0] <- eps; cdf[cdf == 1] <- 1 - eps
  logl <- -Inf; oldlogl <- NaN; v <- 0
  bestlogl <- -Inf; bestf <- f; iters <- 0L

  # MATLAB: while ~(rel<tol) && iters<max && ~isnan(logl). A NaN ratio (first
  # iteration, where oldlogl is -Inf/NaN) counts as "keep going".
  cont <- function() {
    if (is.nan(logl) || iters >= maxIters) return(FALSE)
    rel <- abs(logl - oldlogl) / abs(oldlogl)
    if (!is.finite(rel)) return(TRUE)
    rel >= errorTolerance
  }
  while (cont()) {
    oldlogl <- logl
    pdf <- dnorm(z)
    rc <- pdf / cdf; rw <- pdf / (1 - cdf)
    d <- cpos * (rc * z + rc^2) + cneg * (rw * (-z) + rw^2)
    H <- -crossprod(A, A * d) - regcov            # -A' diag(d) A - regcov
    F <- as.numeric(crossprod(A, cpos * rc - cneg * rw))
    step <- tryCatch(solve(H + diag(diag(H)) * v, F - as.numeric(regcov %*% f)),
                     error = function(e) rep(0, n))
    f <- f - step
    z <- as.numeric(A %*% f)
    cdf <- pnorm(z); cdf[cdf == 0] <- eps; cdf[cdf == 1] <- 1 - eps
    logl <- sum(cpos * log(cdf)) + sum(cneg * log(1 - cdf))
    if (is.finite(logl) && logl > bestlogl) {
      bestlogl <- logl; bestf <- f
    } else {
      v <- max(v * 10, 1e-2)                       # adaptive damping
      f <- bestf; logl <- bestlogl; oldlogl <- logl * 1.1
      z <- as.numeric(A %*% f)
      cdf <- pnorm(z); cdf[cdf == 0] <- eps; cdf[cdf == 1] <- 1 - eps
    }
    iters <- iters + 1L
  }
  bestf
}

# ---- calibrate (port of calibrate.m): mean-variance match to 2k target --------
.paico_calibrate <- function(signal, binAges, target_ages, target_vals,
                             overlap = c(0, 2000)) {
  ov_sig <- which(binAges >= overlap[1] & binAges <= overlap[2])
  if (length(ov_sig) < 3) return(signal)
  part <- signal[ov_sig]
  ov_t <- which(target_ages >= overlap[1] & target_ages <= overlap[2])
  instru <- target_vals[ov_t]
  sp <- stats::sd(part, na.rm = TRUE); si <- stats::sd(instru, na.rm = TRUE)
  if (!is.finite(sp) || sp == 0) return(signal)
  mul <- si / sp
  (signal - mean(part, na.rm = TRUE)) * mul + mean(instru, na.rm = TRUE)
}

run_paico <- function(fts, bandIdx, binvec, binAges, nens,
                      cps_targets = NULL, area_weight, band_weights,
                      apply_reference, cfg = list()) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  N_BANDS <- length(band_weights)
  regcov_scalar <- cfg$paico_reg_param %||% 1

  # PaiCo_12k.m bins at binWidth = 200 yr (far fewer pairwise comparisons than the
  # 100-yr grid). Reconstruct on the 200-yr grid, then interpolate the calibrated
  # band signal onto the shared output grid for the consensus.
  pstep <- 200
  pbinvec <- seq(min(binvec), max(binvec), by = pstep)
  pbinAges <- rowMeans(cbind(pbinvec[-1], pbinvec[-length(pbinvec)]))

  one_member <- function(i) {
    bandMat <- matrix(NA_real_, nrow = length(binAges), ncol = N_BANDS)
    for (b in seq_len(N_BANDS)) {
      sel <- which(bandIdx == b)                       # PaiCo uses all records
      if (length(sel) < 2) next
      binMat <- vapply(fts[sel], function(ts) {
        out <- tryCatch(
          compositeR::sampleEnsembleThenBinTs(ts = ts, binvec = pbinvec, ageVar = "age",
                                              spread = TRUE, alignInterpDirection = FALSE),
          error = function(e) rep(NA_real_, length(pbinAges)))
        if (length(out) != length(pbinAges)) rep(NA_real_, length(pbinAges)) else out
      }, numeric(length(pbinAges)))
      if (is.null(dim(binMat))) binMat <- matrix(binMat, ncol = length(sel))
      binMat[is.nan(binMat)] <- NA

      cc <- .paico_counts(binMat)
      if (is.null(cc)) next
      f <- tryCatch(.paico_optimize(cc$A, cc$cpos, cc$cneg, regcov_scalar),
                    error = function(e) NULL)
      if (is.null(f)) next
      sig <- rep(NA_real_, length(pbinAges))
      present <- which(rowSums(is.finite(binMat)) > 0)   # bins (rows) with any data
      sig[present] <- f[present]
      if (!is.null(cps_targets) && !is.null(cps_targets[[b]])) {
        tgt <- cps_targets[[b]]
        # Calibrate to the target MEDIAN across members (targetMedian.CPS, the
        # series PaiCo_12k.m actually used), not a random member.
        tmed <- apply(tgt$mat, 1, median, na.rm = TRUE)
        sig <- .paico_calibrate(sig, pbinAges, tgt$ages, tmed, overlap = c(0, 2000))
      }
      # interpolate the 200-yr signal onto the shared output grid
      bandMat[, b] <- approx(pbinAges, sig, xout = binAges, rule = 2)$y
    }
    bandMat                                   # full (nbins x N_BANDS) per member
  }

  ncores <- as.integer(cfg$ncores %||% 1)
  if (ncores > 1 && .Platform$OS.type == "unix") {
    cols <- parallel::mclapply(seq_len(nens), one_member, mc.cores = ncores,
                               mc.preschedule = FALSE)
  } else {
    cols <- lapply(seq_len(nens), one_member)
  }
  nb <- length(binAges)
  cols <- lapply(cols, function(m)
    if (is.matrix(m) && all(dim(m) == c(nb, N_BANDS))) m else matrix(NA_real_, nb, N_BANDS))
  rs <- cfg$ref_start %||% 1800; re <- cfg$ref_end %||% 1900
  globalEns <- apply_reference(vapply(cols, area_weight, numeric(nb)), binAges, rs, re)
  bandEns <- lapply(seq_len(N_BANDS), function(b)
    apply_reference(vapply(cols, function(m) m[, b], numeric(nb)), binAges, rs, re))
  list(global = globalEns, bands = bandEns)
}

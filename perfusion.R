# perfusion.R — stress-CMR quantitative perfusion indices (2026-11-13)
#
# Lifted from oht-stress-meta/scripts/stressDB.R (the PhD-thesis definitions, 2024), the
# only copy where MPRi was right. The later oht-hyperemia-harefield copy overwrote MPRi
# with plain MPR and divided stress MBF by the REST RPP; do not use it.
#
#   rest RPP     = rest HR x rest SBP                       (bpm x mmHg)
#   stress RPP   = stress HR x stress SBP
#   rest MBFi    = rest MBF / rest RPP x 1e4                (rest flow corrected for workload)
#   MPR          = stress MBF / rest MBF
#   MPRi         = stress MBF / rest MBFi = MPR x rest RPP / 1e4
#   (alt) MPR / rest RPP x 1e4 is a DIFFERENT quantity, returned as mpr_index_alt only
#   because some older scripts (mpri_splines.R, table4_prognosis.R) used it.
#
# Missing input = NA output (never 0). Physiologically impossible inputs -> NA with a warning.
# Pure function, vectorised; no globals.

perfusion_indices <- function(rest_mbf, stress_mbf, rest_hr = NA_real_, rest_sbp = NA_real_,
                              stress_hr = NA_real_, stress_sbp = NA_real_, mpr = NULL) {
  n <- max(length(rest_mbf), length(stress_mbf), length(rest_hr), length(rest_sbp),
           length(stress_hr), length(stress_sbp))
  rec <- function(v) as.numeric(rep_len(v, n))
  rest_mbf <- rec(rest_mbf); stress_mbf <- rec(stress_mbf)
  rest_hr <- rec(rest_hr); rest_sbp <- rec(rest_sbp)
  stress_hr <- rec(stress_hr); stress_sbp <- rec(stress_sbp)

  guard <- function(v, lo, hi, nm) {
    bad <- !is.na(v) & (v <= lo | v > hi)
    if (any(bad)) warning(sprintf("perfusion_indices: %d implausible %s set to NA", sum(bad), nm))
    v[bad] <- NA_real_; v
  }
  rest_mbf   <- guard(rest_mbf, 0, 10, "rest_mbf")
  stress_mbf <- guard(stress_mbf, 0, 10, "stress_mbf")
  rest_hr    <- guard(rest_hr, 20, 250, "rest_hr")
  stress_hr  <- guard(stress_hr, 20, 250, "stress_hr")
  rest_sbp   <- guard(rest_sbp, 40, 300, "rest_sbp")
  stress_sbp <- guard(stress_sbp, 40, 300, "stress_sbp")

  # If only stress MBF + MPR are known (common in the registry), back-derive rest MBF.
  if (!is.null(mpr)) {
    mpr <- rec(mpr)
    fill <- is.na(rest_mbf) & !is.na(stress_mbf) & !is.na(mpr) & mpr > 0
    rest_mbf[fill] <- stress_mbf[fill] / mpr[fill]
  }

  rest_rpp   <- rest_hr * rest_sbp
  stress_rpp <- stress_hr * stress_sbp
  mpr_calc   <- stress_mbf / rest_mbf
  rest_mbfi  <- rest_mbf / rest_rpp * 1e4
  data.frame(
    rest_rpp = rest_rpp, stress_rpp = stress_rpp,
    mpr = mpr_calc,
    rest_mbfi = rest_mbfi,
    mpri = stress_mbf / rest_mbfi,              # = mpr * rest_rpp / 1e4
    mpr_index_alt = mpr_calc / rest_rpp * 1e4,  # older alternative, sensitivity only
    rpp_ratio = stress_rpp / rest_rpp           # haemodynamic response to the stressor
  )
}

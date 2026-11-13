# cmr_plots.R — CMR figure helpers (2026-11-13)
#
#   aha_bullseye_polygons()  AHA-16/17 segment polygons (data frame for ggplot)
#   plot_bullseye()          per-segment values -> bullseye
#   plot_spline_hr()         HR curve (95% CI) for one continuous term of a coxph fit
#
# Replaces:
#   * aha-16-segment (AHA16::bullseye) and its acm-harefield / myosport forks — base
#     graphics, needed {fields}, values below the scale gave a blank segment, and the
#     forks indexed a 264-colour palette with a 255/256-colour one so the top ~3 % of
#     the scale also came out blank. Geometry here is the rbhh-cmr3 P13 ring geometry
#     (R/projects/p13_mbf_normals/03_segment_reference.R).
#   * plotHR() (farm-splines / hcm-lge / oht-stress-meta / myocarditis-bristol; originally
#     R. Seifert & M. Gordon, later Greg::plotHR) — re-read the data via eval(model$call$data),
#     mis-handled interaction terms, changed par() globally. This version takes `data`
#     explicitly, holds the other covariates at stated values and returns a ggplot.
#
# Needs: ggplot2 (+ survival for plot_spline_hr).

# ---- bullseye ---------------------------------------------------------------

#' AHA segment polygons. Orientation: anterior at the top, septum on the LEFT
#' (standard short-axis display). Basal ring outer (1-6, 60 deg each, seg 1 at top),
#' mid (7-12), apical (13-16, 90 deg each), optional apex 17 as the centre disc.
aha_bullseye_polygons <- function(apex17 = FALSE, n = 30) {
  ring <- function(segs, r0, r1, start, width) do.call(rbind, lapply(seq_along(segs), function(i) {
    a <- seq(start + (i - 1) * width, start + i * width, length.out = n) * pi / 180
    data.frame(segment = segs[i],
               x = c(r1 * cos(a), rev(r0 * cos(a))), y = c(r1 * sin(a), rev(r0 * sin(a))),
               lx = (r0 + r1) / 2 * cos(mean(a)), ly = (r0 + r1) / 2 * sin(mean(a)))
  }))
  r_in <- if (apex17) .35 else 0
  out <- rbind(ring(1:6, 2, 3, 60, 60), ring(7:12, 1, 2, 60, 60), ring(13:16, r_in, 1, 45, 90))
  if (apex17) {
    a <- seq(0, 2 * pi, length.out = 4 * n)
    out <- rbind(out, data.frame(segment = 17, x = r_in * cos(a), y = r_in * sin(a), lx = 0, ly = 0))
  }
  out
}

#' Bullseye of per-segment values.
#' @param values numeric, either length 16/17 in segment order, or named by segment
#'   number ("1".."17"). NA = segment not known -> drawn in `na_fill`, never as 0.
#' @param limits colour-scale limits; values outside are SQUISHED to the ends (not blank).
#' @param palette vector of colours for a continuous gradient (low -> high).
#' @param labels TRUE to print the value in each segment; or a character vector.
plot_bullseye <- function(values, limits = NULL, palette = c("#e6f1fb", "#185fa5"),
                          title = NULL, fill_name = NULL, labels = TRUE, digits = 2,
                          na_fill = "grey85", territories = FALSE, base_size = 10) {
  if (is.null(names(values))) {
    if (!length(values) %in% c(16, 17)) stop("plot_bullseye: give 16/17 values or name them by segment")
    names(values) <- seq_along(values)
  }
  segs <- as.integer(names(values))
  if (any(is.na(segs)) || any(!segs %in% 1:17)) stop("plot_bullseye: names must be segment numbers 1-17")
  poly <- aha_bullseye_polygons(apex17 = 17 %in% segs)
  poly$value <- unname(values[as.character(poly$segment)])
  if (is.null(limits)) limits <- range(values, na.rm = TRUE)
  lab <- unique(poly[c("segment", "lx", "ly", "value")])
  if (isTRUE(labels)) {
    lab$txt <- ifelse(is.na(lab$value), "", formatC(lab$value, format = "f", digits = digits))
  } else if (is.character(labels)) {
    lab$txt <- labels[match(lab$segment, seq_along(labels))]
  } else lab <- lab[0, ]

  p <- ggplot2::ggplot(poly, ggplot2::aes(x, y)) +
    ggplot2::geom_polygon(ggplot2::aes(group = segment, fill = value), colour = "white", linewidth = .6) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(colours = palette, limits = limits, oob = scales::squish,
                                  na.value = na_fill, name = fill_name) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::labs(title = title)
  if (nrow(lab)) p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(lx, ly, label = txt),
                                             size = base_size / ggplot2::.pt * .8)
  if (territories) {
    # typical dominance: LAD anterior/septal (top), RCA inferior (bottom), LCX lateral (right)
    p <- p + ggplot2::annotate("text", x = c(0, 0, 3.55), y = c(3.35, -3.35, 0),
                               label = c("LAD", "RCA", "LCX"), size = base_size / ggplot2::.pt,
                               fontface = "bold", colour = "grey30")
  }
  p
}

# ---- spline hazard ratio ----------------------------------------------------

#' HR curve for one continuous variable of a coxph fit, relative to `ref`.
#' Works with linear terms, splines::ns(), survival::pspline() and other basis terms
#' that keep their knots in the model terms. Other covariates are held at `at`
#' (default: median of numerics, most frequent level of factors/characters, taken
#' from `data`) — they cancel in the HR unless they interact with `var`.
#' @param data the data the model was fitted on (passed explicitly; never re-evaluated).
#' @return a ggplot; the curve data are in attr(p, "curve").
plot_spline_hr <- function(fit, data, var, ref = NULL, at = list(), n = 200,
                           probs = c(.025, .975), conf = 0.95, log_y = TRUE,
                           rug = TRUE, xlab = var, ylab = "Hazard ratio",
                           colour = "#005694", fill = "#DEEBF7") {
  if (!inherits(fit, "coxph")) stop("plot_spline_hr: needs a coxph fit")
  tt <- stats::delete.response(stats::terms(fit))
  if (any(grepl("strata\\(", attr(tt, "term.labels")))) stop("plot_spline_hr: strata() terms not supported")
  vars <- all.vars(tt)
  if (!var %in% vars) stop("plot_spline_hr: '", var, "' is not in the model")
  x <- data[[var]]; x <- x[is.finite(x)]
  if (is.null(ref)) ref <- stats::median(x)
  grid_x <- seq(stats::quantile(x, probs[1]), stats::quantile(x, probs[2]), length.out = n)

  typical <- function(v) {
    if (is.numeric(v)) return(stats::median(v, na.rm = TRUE))
    tab <- table(v); lv <- names(tab)[which.max(tab)]
    if (is.factor(v)) factor(lv, levels = levels(v)) else lv
  }
  base <- lapply(setdiff(vars, var), function(v) if (!is.null(at[[v]])) at[[v]] else typical(data[[v]]))
  names(base) <- setdiff(vars, var)
  mk <- function(xx) { d <- as.data.frame(base, stringsAsFactors = FALSE)[rep(1, length(xx)), , drop = FALSE]
                       d[[var]] <- xx; d }
  X <- function(d) { mf <- stats::model.frame(tt, d, xlev = fit$xlevels)
                     m <- stats::model.matrix(tt, mf); m[, colnames(m) != "(Intercept)", drop = FALSE] }
  Xg <- X(mk(grid_x)); Xr <- X(mk(ref))
  b <- stats::coef(fit); V <- stats::vcov(fit)
  if (ncol(Xg) != length(b)) stop("plot_spline_hr: design matrix does not match the coefficients")
  keep <- !is.na(b); D <- sweep(Xg[, keep, drop = FALSE], 2, Xr[1, keep])
  lp <- drop(D %*% b[keep]); se <- sqrt(pmax(0, rowSums((D %*% V[keep, keep]) * D)))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  curve <- data.frame(x = grid_x, hr = exp(lp), lo = exp(lp - z * se), hi = exp(lp + z * se))

  p <- ggplot2::ggplot(curve, ggplot2::aes(x, hr)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = fill) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    ggplot2::geom_line(colour = colour, linewidth = 1) +
    ggplot2::labs(x = xlab, y = ylab,
                  caption = sprintf("Reference %s = %s; %g%% CI", var, format(signif(ref, 3)), 100 * conf)) +
    ggplot2::theme_classic()
  if (log_y) p <- p + ggplot2::scale_y_log10()
  if (rug) p <- p + ggplot2::geom_rug(data = data.frame(x = x[x >= min(grid_x) & x <= max(grid_x)]),
                                      ggplot2::aes(x = x), inherit.aes = FALSE, alpha = .15, length = ggplot2::unit(.02, "npc"))
  attr(p, "curve") <- curve
  p
}

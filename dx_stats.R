# dx_stats.R — diagnostic-test helpers: likelihood ratios, post-test probability,
# Fagan nomogram (2026-11-13)
#
# Replaces jacc-cmd-letter/leafplot.R (dis.pos.given.test.pos/neg, inputs in %) and
# nomogrammer.R (A. Chekroud's MIT-licensed nomogrammer, vendored there without
# attribution; it also called theme_set() and so restyled every later plot).
# Everything here takes PROBABILITIES on the 0-1 scale, never percent.
#
# Needs: ggplot2 (fagan_nomogram only).

.check_prob <- function(p, nm, open = TRUE) {
  if (any(!is.na(p) & (p < 0 | p > 1))) stop(nm, " must be a probability in [0, 1] (did you pass a %?)")
  if (open && any(!is.na(p) & (p == 0 | p == 1))) stop(nm, " must be strictly between 0 and 1")
  invisible(p)
}

#' Positive and negative likelihood ratios from sensitivity and specificity.
likelihood_ratios <- function(sens, spec) {
  .check_prob(sens, "sens", open = FALSE); .check_prob(spec, "spec", open = FALSE)
  data.frame(plr = sens / (1 - spec), nlr = (1 - sens) / spec)
}

#' Post-test probability from pre-test probability and a likelihood ratio
#' (Bayes in odds form). Vectorised.
post_test_prob <- function(pre, lr) {
  .check_prob(pre, "pre")
  if (any(!is.na(lr) & lr < 0)) stop("lr must be >= 0")
  o <- pre / (1 - pre) * lr
  o / (1 + o)
}

#' Post-test probabilities after a positive and a negative test.
post_test <- function(pre, sens, spec) {
  lr <- likelihood_ratios(sens, spec)
  data.frame(pre = pre, plr = lr$plr, nlr = lr$nlr,
             post_pos = post_test_prob(pre, lr$plr),
             post_neg = post_test_prob(pre, lr$nlr))
}

#' Fagan nomogram (ggplot). Give sens + spec, or plr + nlr.
#' Geometry: left axis = -logit(pre), right = logit(post), middle = log(LR) / 2, so each
#' straight line through the LR mark lands on the correct post-test probability.
fagan_nomogram <- function(pre, sens = NULL, spec = NULL, plr = NULL, nlr = NULL,
                           colours = c(positive = "#ad1122", negative = "#005694"),
                           show_values = TRUE, base_size = 11) {
  .check_prob(pre, "pre")
  if (length(pre) != 1) stop("fagan_nomogram: one pre-test probability")
  if (!is.null(sens) && !is.null(spec)) {
    lr <- likelihood_ratios(sens, spec); plr <- lr$plr; nlr <- lr$nlr
  } else if (is.null(plr) || is.null(nlr)) {
    stop("fagan_nomogram: give sens + spec, or plr + nlr")
  }
  lg <- function(p) log10(p / (1 - p))
  probs <- c(.001, .005, .01, .02, .05, .1, .2, .3, .5, .7, .8, .9, .95, .98, .99, .995, .999)
  lrs <- c(.001, .01, .1, .2, .5, 1, 2, 5, 10, 100, 1000)
  lim <- range(lg(probs))
  pos <- post_test_prob(pre, plr); neg <- post_test_prob(pre, nlr)
  lines <- data.frame(
    test = rep(c("positive", "negative"), each = 2), x = rep(c(0, 1), 2),
    y = c(-lg(pre), lg(pos), -lg(pre), lg(neg)))
  fmtp <- function(p) paste0(trimws(formatC(100 * p, format = "fg", digits = 3)), "%")
  ticks <- function(x, y, lab, hjust) data.frame(x = x, y = y, lab = lab, hjust = hjust)
  tk <- rbind(
    ticks(0, -lg(probs), fmtp(probs), 1.15),
    ticks(1, lg(probs), fmtp(probs), -0.15),
    ticks(0.5, log10(lrs) / 2, format(lrs, drop0trailing = TRUE, scientific = FALSE, trim = TRUE), -0.25))
  tk <- tk[tk$y >= lim[1] - 1e-9 & tk$y <= lim[2] + 1e-9, ]

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = c(0, .5, 1), xend = c(0, .5, 1), y = lim[1], yend = lim[2]),
                          linewidth = .4) +
    ggplot2::geom_segment(data = tk, ggplot2::aes(x = x - .012, xend = x + .012, y = y, yend = y), linewidth = .3) +
    ggplot2::geom_text(data = tk, ggplot2::aes(x = x, y = y, label = lab, hjust = hjust),
                       size = base_size / ggplot2::.pt * .8) +
    ggplot2::geom_line(data = lines, ggplot2::aes(x = x, y = y, colour = test), linewidth = 1) +
    ggplot2::annotate("text", x = c(0, .5, 1), y = lim[2] + .35,
                      label = c("Pre-test\nprobability", "Likelihood\nratio", "Post-test\nprobability"),
                      size = base_size / ggplot2::.pt, fontface = "bold", lineheight = .9) +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::scale_y_continuous(limits = c(lim[1], lim[2] + .6)) +
    ggplot2::scale_x_continuous(limits = c(-.2, 1.2)) +
    ggplot2::theme_void(base_size = base_size)
  if (show_values) {
    p <- p + ggplot2::labs(caption = sprintf(
      "Pre-test %s. LR+ %.2f -> post-test %s; LR- %.2f -> post-test %s",
      fmtp(pre), plr, fmtp(pos), nlr, fmtp(neg)))
  }
  p
}

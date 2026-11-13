# Tests for perfusion.R, dx_stats.R, cmr_plots.R.
# Run from the repo root:  Rscript -e 'testthat::test_dir("tests")'
library(testthat)
root <- if (file.exists("perfusion.R")) "." else ".."
for (f in c("perfusion.R", "dx_stats.R", "cmr_plots.R")) source(file.path(root, f))

test_that("perfusion_indices: thesis definitions", {
  r <- perfusion_indices(rest_mbf = 1, stress_mbf = 3, rest_hr = 60, rest_sbp = 120,
                         stress_hr = 90, stress_sbp = 130)
  expect_equal(r$rest_rpp, 7200); expect_equal(r$stress_rpp, 11700)
  expect_equal(r$mpr, 3)
  expect_equal(r$rest_mbfi, 1 / 7200 * 1e4)
  expect_equal(r$mpri, 3 * 7200 / 1e4)                 # = MPR x RPP / 1e4 (P18 lib_p18.R)
  expect_equal(r$mpr_index_alt, 3 / 7200 * 1e4)
})

test_that("perfusion_indices: missing stays NA, implausible -> NA, rest back-derived from MPR", {
  r <- perfusion_indices(rest_mbf = c(1, NA), stress_mbf = c(2, 2.4), rest_hr = c(NA, 60), rest_sbp = 120)
  expect_equal(r$mpr, c(2, NA)); expect_true(all(is.na(r$mpri)))
  expect_warning(r2 <- perfusion_indices(rest_mbf = 0, stress_mbf = 2), "implausible")
  expect_true(is.na(r2$mpr))
  r3 <- perfusion_indices(rest_mbf = NA, stress_mbf = 2.4, mpr = 2, rest_hr = 60, rest_sbp = 100)
  expect_equal(r3$mpr, 2); expect_equal(r3$mpri, 2 * 6000 / 1e4)
})

test_that("likelihood ratios and post-test probability", {
  lr <- likelihood_ratios(0.8, 0.9)
  expect_equal(lr$plr, 8); expect_equal(lr$nlr, 0.2 / 0.9)
  expect_equal(post_test_prob(0.5, 8), 8 / 9)
  expect_equal(post_test_prob(0.2, 1), 0.2)
  pt <- post_test(0.3, 0.8, 0.9)
  expect_equal(pt$post_pos, post_test_prob(0.3, 8))
  expect_error(post_test_prob(30, 2), "probability")
  old_pos <- function(pre, sens, spec) {               # jacc-cmd-letter/leafplot.R (inputs in %)
    (sens * pre) / (sens * pre + (100 - spec) * (100 - pre)) * 100
  }
  expect_equal(100 * post_test(0.25, 0.7, 0.85)$post_pos, old_pos(25, 70, 85))
})

test_that("fagan_nomogram builds; lines land on the post-test probabilities", {
  p <- fagan_nomogram(0.3, sens = 0.8, spec = 0.9)
  expect_s3_class(p, "ggplot")
  ld <- ggplot2::layer_data(p, 4)                       # the two test lines
  lg <- function(q) log10(q / (1 - q))
  right <- sort(ld$y[ld$x == 1])
  expect_equal(right, sort(lg(c(post_test_prob(0.3, 8), post_test_prob(0.3, 0.2 / 0.9)))))
  mid <- (ld$y[ld$x == 0][1] + ld$y[ld$x == 1][1]) / 2   # straight line crosses the middle at log10(LR)/2
  expect_true(any(abs(mid - log10(c(8, 0.2 / 0.9)) / 2) < 1e-12))
  expect_error(fagan_nomogram(0.3), "sens")
})

test_that("bullseye geometry: orientation and segment count", {
  g <- aha_bullseye_polygons()
  expect_setequal(unique(g$segment), 1:16)
  c1 <- unique(g[g$segment == 1, c("lx", "ly")]); c2 <- unique(g[g$segment == 2, c("lx", "ly")])
  c4 <- unique(g[g$segment == 4, c("lx", "ly")]); c6 <- unique(g[g$segment == 6, c("lx", "ly")])
  expect_true(abs(c1$lx) < 1e-9 && c1$ly > 0)           # 1 basal anterior: top
  expect_true(c2$lx < 0)                                # 2 anteroseptal: left
  expect_true(abs(c4$lx) < 1e-9 && c4$ly < 0)           # 4 inferior: bottom
  expect_true(c6$lx > 0)                                # 6 anterolateral: right
  expect_setequal(unique(aha_bullseye_polygons(apex17 = TRUE)$segment), 1:17)
})

test_that("plot_bullseye: NA grey, out-of-range squished (not blank)", {
  v <- c(seq(0.5, 4, length.out = 15), NA)
  p <- plot_bullseye(v, limits = c(1, 3))
  ld <- ggplot2::layer_data(p, 1)
  fills <- tapply(ld$fill, ld$group, `[`, 1)
  expect_equal(sum(fills == "grey85"), 1)               # only the NA segment
  expect_false(any(is.na(fills)))
  seg_fill <- function(s) ld$fill[ld$group == s][1]     # groups follow segment order
  expect_equal(toupper(seg_fill(1)), "#E6F1FB")         # 0.5 below the limit -> low-end colour
  expect_equal(toupper(seg_fill(15)), "#185FA5")        # 4.0 above the limit -> high-end colour
  expect_error(plot_bullseye(1:5), "16/17")
  expect_s3_class(plot_bullseye(c(`3` = 2, `17` = 1)), "ggplot")
})

test_that("plot_spline_hr: linear term reproduces exp(beta * dx); spline fits run", {
  suppressMessages({library(survival); library(splines)})
  d <- na.omit(lung[, c("time", "status", "age", "sex", "ph.ecog")])
  d$ecog <- factor(d$ph.ecog)
  f <- coxph(Surv(time, status) ~ age + sex + ecog, data = d)
  p <- plot_spline_hr(f, d, "age", ref = 60)
  cv <- attr(p, "curve")
  expect_equal(cv$hr, exp(coef(f)[["age"]] * (cv$x - 60)), tolerance = 1e-10)
  se <- sqrt(vcov(f)["age", "age"]) * abs(cv$x - 60)
  expect_equal(cv$hi, exp(coef(f)[["age"]] * (cv$x - 60) + qnorm(.975) * se), tolerance = 1e-8)
  f2 <- coxph(Surv(time, status) ~ ns(age, 3) + sex, data = d)
  cv2 <- attr(plot_spline_hr(f2, d, "age"), "curve")
  expect_equal(cv2$hr[which.min(abs(cv2$x - median(d$age)))], 1, tolerance = 0.05)
  f3 <- coxph(Surv(time, status) ~ pspline(age, df = 3) + sex, data = d)
  expect_s3_class(plot_spline_hr(f3, d, "age"), "ggplot")
  expect_error(plot_spline_hr(f, d, "mpr"), "not in the model")
})

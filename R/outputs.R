## Tidy outputs and base-R plots for manuscripts.
## All base R: no ggplot2, flextable, or other dependencies. Each function is
## short enough to read and audit in one sitting. Returns objects (data frames,
## or draws to the active device); never writes files.

#' @importFrom graphics abline lines par points rect
#' @importFrom grDevices adjustcolor
NULL

# ---- Predictive: tidy data frame -------------------------------------------

#' Coerce an `hla_ci` result to a data frame
#'
#' Returns the metric / estimate / lower / upper rows as a plain data frame,
#' ready to pipe into any table formatter (`flextable`, `gt`, `knitr::kable`,
#' `huxtable`) in a manuscript. The package deliberately does not depend on a
#' table package; it produces the numbers and leaves formatting to the user.
#'
#' @param x An object from [ppv_npv_nnt_ci()].
#' @param ... Ignored.
#' @return A data frame with columns `metric`, `estimate`, `lower`, `upper`.
#' @examples
#' r <- ppv_npv_nnt_ci(12, 29, 2, 49, incidence = c(1/2000, 1/500),
#'                     population_carriers = 3255, population_n = 54249, seed = 1)
#' as.data.frame(r)
#' @export
as.data.frame.hla_ci <- function(x, ...) {
  data.frame(
    metric   = x$ci$metric,
    estimate = x$ci$estimate,
    lower    = x$ci$lower,
    upper    = x$ci$upper,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# ---- Predictive: signature plot --------------------------------------------

#' Plot carrier risk (or NNT) versus incidence
#'
#' The figure that makes the package's central point: for a fixed allele the
#' clinical utility is driven by the phenotype incidence, not the allele. Plots
#' carrier risk (PPV) or NNT against incidence on log axes, with the typical
#' SCAR incidence range shaded, and marks the current estimate.
#'
#' Pure base-R graphics; draws to the active device and returns `NULL`
#' invisibly.
#'
#' @param sensitivity Proportion of cases carrying the allele.
#' @param carriage Population allele carriage frequency.
#' @param metric Either `"carrier_risk"` (PPV) or `"nnt"`.
#' @param incidence_lo,incidence_hi Incidence axis range (as proportions).
#' @param mark_incidence Optional incidence to mark with a point.
#' @param scar_band Length-2 incidence range to shade (default 1e-4 to 1e-3).
#' @param col Line/point colour.
#' @return Invisibly `NULL`; called for its plot.
#' @examples
#' plot_carrier_risk_vs_incidence(sensitivity = 0.414, carriage = 0.06,
#'                                mark_incidence = 1/1000)
#' @export
plot_carrier_risk_vs_incidence <- function(sensitivity, carriage,
                                           metric = c("carrier_risk", "nnt"),
                                           incidence_lo = 1e-5, incidence_hi = 1e-1,
                                           mark_incidence = NULL,
                                           scar_band = c(1e-4, 1e-3),
                                           col = "#3F6699") {
  metric <- match.arg(metric)
  I <- 10 ^ seq(log10(incidence_lo), log10(incidence_hi), length.out = 200)
  if (metric == "carrier_risk") {
    y <- pmin(1, sensitivity * I / carriage); ylab <- "Carrier risk (PPV)"; logaxes <- "x"
  } else {
    y <- 1 / (sensitivity * I); ylab <- "Number needed to test"; logaxes <- "xy"
  }
  plot(I, y, type = "n", log = logaxes, xlab = "Incidence among exposed", ylab = ylab,
       main = sprintf("%s vs incidence (sensitivity %.0f%%, carriage %.1f%%)",
                      ylab, 100 * sensitivity, 100 * carriage))
  if (length(scar_band) == 2) {
    u <- par("usr"); ylo <- if (logaxes == "xy") 10^u[3] else u[3]
    yhi <- if (logaxes == "xy") 10^u[4] else u[4]
    rect(scar_band[1], ylo, scar_band[2], yhi, col = adjustcolor(col, 0.08), border = NA)
  }
  lines(I, y, lwd = 2.5, col = col)
  if (!is.null(mark_incidence)) {
    ym <- if (metric == "carrier_risk") min(1, sensitivity * mark_incidence / carriage)
          else 1 / (sensitivity * mark_incidence)
    points(mark_incidence, ym, pch = 19, col = col, cex = 1.3)
  }
  invisible(NULL)
}

# ---- Power: scenario grid and curve ----------------------------------------

#' Minimum detectable OR across a grid of designs
#'
#' Tabulates the minimum detectable odds ratio (and power at a reference OR)
#' over combinations of case count, control:case ratio, carriage, and alpha.
#' Returns a tidy data frame for a power-justification table or supplement.
#'
#' @param n_case Vector of case counts.
#' @param carriage Vector of control/population carriage frequencies.
#' @param ratio Vector of controls-per-case.
#' @param alpha Vector of significance thresholds.
#' @param power Target power for the detectable OR (default 0.80).
#' @param ref_or Reference OR at which to also report power (default 10).
#' @return A data frame with one row per scenario: `n_case`, `ratio`, `n_ctrl`,
#'   `carriage`, `alpha`, `mdor`, and `power_at_ref_or`.
#' @examples
#' power_grid(n_case = c(20, 29, 50), carriage = 0.04, ratio = 2, alpha = 0.05)
#' @export
power_grid <- function(n_case, carriage, ratio = 1, alpha = 0.05,
                       power = 0.80, ref_or = 10) {
  g <- expand.grid(n_case = n_case, carriage = carriage, ratio = ratio,
                   alpha = alpha, KEEP.OUT.ATTRS = FALSE)
  g$n_ctrl <- round(g$n_case * g$ratio)
  g$mdor <- mapply(function(nc, f, a, nt)
    min_detectable_or(nc, f, alpha = a, n_ctrl = nt, power = power),
    g$n_case, g$carriage, g$alpha, g$n_ctrl)
  g$power_at_ref_or <- mapply(function(nc, f, a, nt)
    fisher_power(nc, f, or = ref_or, alpha = a, n_ctrl = nt),
    g$n_case, g$carriage, g$alpha, g$n_ctrl)
  g[, c("n_case", "ratio", "n_ctrl", "carriage", "alpha", "mdor", "power_at_ref_or")]
}

#' Plot the exact power curve versus odds ratio
#'
#' Power as a function of the odds ratio at a fixed sample size, with the target
#' power line and the minimum detectable OR marked. The figure for the
#' STROBE/STREGA minimum-detectable-effect statement. Pure base-R graphics.
#'
#' @inheritParams fisher_power
#' @param target Target power line to draw (default 0.80).
#' @param or_hi Upper OR for the x-axis (default 40).
#' @param mark_or Optional observed OR to mark.
#' @param col Line/marker colour.
#' @return Invisibly the minimum detectable OR at `target`.
#' @examples
#' plot_power_curve(n_case = 29, carriage = 0.041, ratio = 2)
#' @export
plot_power_curve <- function(n_case, carriage, alpha = 0.05, ratio = 1,
                             n_ctrl = NULL, target = 0.80, or_hi = 40,
                             mark_or = NULL, col = "#3F6699") {
  if (is.null(n_ctrl)) n_ctrl <- round(n_case * ratio)
  reject <- build_reject(n_case, n_ctrl, alpha)
  ors <- seq(1, or_hi, by = 0.5)
  pw <- vapply(ors, function(o)
    fisher_power(n_case, carriage, o, alpha, n_ctrl = n_ctrl, reject = reject), numeric(1))
  plot(ors, pw, type = "l", lwd = 2.5, col = col, ylim = c(0, 1),
       xlab = "Odds ratio", ylab = "Power",
       main = sprintf("Exact power (n_case = %d, n_ctrl = %d, alpha = %.3g)",
                      n_case, n_ctrl, alpha))
  abline(h = target, lty = 2, col = "grey50")
  md <- min_detectable_or(n_case, carriage, alpha, n_ctrl = n_ctrl, power = target, or_max = or_hi)
  if (!is.na(md)) abline(v = md, lty = 3, col = "grey50")
  if (!is.null(mark_or)) {
    points(mark_or, fisher_power(n_case, carriage, mark_or, alpha, n_ctrl = n_ctrl, reject = reject),
           pch = 19, col = col, cex = 1.3)
  }
  invisible(md)
}

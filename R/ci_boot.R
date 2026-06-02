## Internal helpers -----------------------------------------------------------

# Draw incidence samples. `incidence` may be:
#   - a scalar: treated as fixed (no incidence uncertainty)
#   - length-2 c(lo, hi): log-uniform between lo and hi
#   - a function(n): user-supplied sampler returning n draws
.draw_incidence <- function(incidence, n) {
  if (is.function(incidence)) return(incidence(n))
  if (length(incidence) == 1L) return(rep(incidence, n))
  if (length(incidence) == 2L) {
    return(exp(stats::runif(n, log(incidence[1]), log(incidence[2]))))
  }
  stop("`incidence` must be a scalar, a length-2 range c(lo, hi), or a function(n).",
       call. = FALSE)
}

# Point incidence for the central estimate (geometric mean of a range).
.incidence_point <- function(incidence) {
  if (is.function(incidence)) return(stats::median(incidence(10000L)))
  if (length(incidence) == 2L) return(exp(mean(log(incidence))))
  incidence
}

# Exact Clopper-Pearson interval for x successes in n (base R).
.clopper_pearson <- function(x, n, conf.level) {
  as.numeric(stats::binom.test(x, n, conf.level = conf.level)$conf.int)
}

#' Confidence intervals for PPV, NPV, NNT, and carrier risk
#'
#' Propagates the sampling uncertainty in sensitivity and specificity, and the
#' (usually dominant) uncertainty in incidence, into interval estimates for the
#' derived clinical-utility metrics. Sensitivity and specificity are drawn from
#' Jeffreys Beta distributions (`Beta(x + 0.5, n - x + 0.5)`), a smoothed
#' parametric bootstrap that, unlike resampling raw counts, stays well-behaved
#' when a cell is zero (e.g. zero carriers among matched controls). Incidence is
#' an explicit uncertain input, because for rare reactions it dominates the
#' uncertainty in PPV and NNT.
#'
#' PPV/NPV use the specificity-based form
#' \eqn{PPV = (Se\,I)/(Se\,I + (1-Sp)(1-I))} (matching [ppv_specificity()] and
#' [predictive_ci()]). For the carriage-based PPV used by [ppv_npv_nnt()],
#' supply `population_carriage`; its CI is then added too.
#'
#' Point estimates are evaluated at the maximum-likelihood inputs (not the
#' bootstrap mean) to avoid resampling bias. NNT intervals use the percentile
#' method because NNT is a skewed ratio for which a symmetric interval would
#' misrepresent the uncertainty. Sensitivity and specificity are reported with
#' exact Clopper-Pearson intervals rather than bootstrap percentiles.
#'
#' @param cases_carrier,cases_n Carriers among cases, and total cases.
#' @param controls_carrier,controls_n Carriers among controls, and total controls.
#' @param incidence Incidence of the reaction among exposed: a scalar (fixed),
#'   a length-2 range `c(lo, hi)` (log-uniform), or a `function(n)` sampler.
#' @param population_carriers,population_n Population carriage as a count and its
#'   sample size (e.g. carriers and total genotyped in a biobank). Supplying both
#'   adds the carriage-based PPV with an interval that propagates carriage
#'   uncertainty. This is the preferred way to specify carriage.
#' @param population_carriage Alternative scalar carriage (e.g. `0.06`) for when
#'   only a point estimate is available with no usable denominator; held fixed.
#'   Mutually exclusive with `population_carriers`/`population_n`.
#' @param n_boot Number of Monte Carlo replicates (default 10000).
#' @param conf.level Confidence level (default 0.95).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class `hla_ci`: a list with `estimate` (named point
#'   estimates), `ci` (a data frame of lower/upper bounds per metric),
#'   `conf.level`, `n_boot`, `seed`, `method`, and any `warnings`.
#'
#' @examples
#' # Konvinse vancomycin DRESS: 19/23 cases, 0/46 controls; incidence uncertain
#' ppv_npv_nnt_ci(cases_carrier = 19, cases_n = 23,
#'                controls_carrier = 0, controls_n = 46,
#'                incidence = c(1/100, 1/40),
#'                population_carriers = 3418, population_n = 54249, seed = 1)
#' @export
ppv_npv_nnt_ci <- function(cases_carrier, cases_n,
                           controls_carrier, controls_n,
                           incidence,
                           population_carriers = NULL, population_n = NULL,
                           population_carriage = NULL,
                           n_boot = 10000L, conf.level = 0.95,
                           seed = NULL) {
  stopifnot(cases_carrier >= 0, cases_carrier <= cases_n, cases_n >= 1,
            controls_carrier >= 0, controls_carrier <= controls_n, controls_n >= 1,
            n_boot >= 100, conf.level > 0, conf.level < 1)

  # Resolve carriage specification: count+n (preferred), scalar, or none.
  have_counts <- !is.null(population_carriers) || !is.null(population_n)
  have_scalar <- !is.null(population_carriage)
  if (have_counts && have_scalar) {
    stop("Specify carriage either as population_carriers + population_n, or as ",
         "the scalar population_carriage, not both.", call. = FALSE)
  }
  if (have_counts && (is.null(population_carriers) || is.null(population_n))) {
    stop("Supply both population_carriers and population_n (or use the scalar ",
         "population_carriage instead).", call. = FALSE)
  }
  if (have_counts) {
    stopifnot(population_carriers >= 0, population_carriers <= population_n,
              population_n >= 1)
  }

  if (!is.null(seed)) set.seed(as.integer(seed))

  a   <- 1 - conf.level
  qlo <- a / 2; qhi <- 1 - a / 2
  pct <- function(x) stats::quantile(x, c(qlo, qhi), names = FALSE, na.rm = TRUE)

  se_hat <- cases_carrier / cases_n
  sp_hat <- 1 - controls_carrier / controls_n
  I_pt   <- .incidence_point(incidence)

  ppv_f <- function(se, sp, I) (se * I) / (se * I + (1 - sp) * (1 - I))
  npv_f <- function(se, sp, I) (sp * (1 - I)) / ((1 - se) * I + sp * (1 - I))
  nnt_f <- function(se, I) 1 / (se * I)

  # Smoothed bootstrap draws (Jeffreys Beta) + incidence sampler
  se_b <- stats::rbeta(n_boot, cases_carrier + 0.5, (cases_n - cases_carrier) + 0.5)
  sp_b <- stats::rbeta(n_boot, (controls_n - controls_carrier) + 0.5, controls_carrier + 0.5)
  I_b  <- .draw_incidence(incidence, n_boot)

  ppv_b <- ppv_f(se_b, sp_b, I_b)
  npv_b <- npv_f(se_b, sp_b, I_b)
  nnt_b <- nnt_f(se_b, I_b)

  estimate <- c(
    sensitivity        = se_hat,
    specificity        = sp_hat,
    ppv_study_controls = ppv_f(se_hat, sp_hat, I_pt),
    npv                = npv_f(se_hat, sp_hat, I_pt),
    nnt                = nnt_f(se_hat, I_pt)
  )

  ci <- rbind(
    sensitivity        = .clopper_pearson(cases_carrier, cases_n, conf.level),
    specificity        = .clopper_pearson(controls_n - controls_carrier, controls_n, conf.level),
    ppv_study_controls = pct(ppv_b),
    npv                = pct(npv_b),
    nnt                = pct(nnt_b)
  )

  # Optional carriage-based PPV (the ppv_npv_nnt() form)
  if (have_counts || have_scalar) {
    if (have_counts) {
      f_b  <- stats::rbeta(n_boot, population_carriers + 0.5,
                           (population_n - population_carriers) + 0.5)
      f_pt <- population_carriers / population_n
    } else {
      f_b  <- rep(population_carriage, n_boot)
      f_pt <- population_carriage
    }
    ppv_carr_b <- pmin(1, se_b * I_b / f_b)
    estimate["ppv_population"] <- min(1, se_hat * I_pt / f_pt)
    ci <- rbind(ci, ppv_population = pct(ppv_carr_b))
  }

  # Guards / warnings
  warns <- character(0)
  min_cell <- min(cases_carrier, cases_n - cases_carrier,
                  controls_carrier, controls_n - controls_carrier)
  if (min_cell < 5) {
    warns <- c(warns, sprintf(
      "Sparse data (smallest cell = %d): intervals are prior-regularized (Jeffreys) and the point estimate for the zero/near-zero rate is not pinned by the data.",
      min_cell))
  }
  if (!is.function(incidence) && length(incidence) == 1L) {
    warns <- c(warns,
      "Incidence treated as fixed: PPV/NPV/NNT intervals reflect only sampling error in sensitivity/specificity and understate total uncertainty. Pass incidence = c(lo, hi) to propagate incidence uncertainty.")
  }
  for (w in warns) warning(w, call. = FALSE)

  ci_df <- data.frame(
    metric = rownames(ci),
    estimate = estimate[rownames(ci)],
    lower = ci[, 1],
    upper = ci[, 2],
    row.names = NULL
  )

  structure(list(
    estimate   = estimate,
    ci         = ci_df,
    conf.level = conf.level,
    n_boot     = n_boot,
    seed       = seed,
    method     = "Jeffreys-Beta smoothed bootstrap; exact Clopper-Pearson for Se/Sp; incidence per `incidence`",
    warnings   = warns
  ), class = "hla_ci")
}

#' @export
print.hla_ci <- function(x, ...) {
  cat(sprintf("PPV/NPV/NNT with %.0f%% intervals (%d replicates%s)\n",
              100 * x$conf.level, x$n_boot,
              if (!is.null(x$seed)) sprintf(", seed %d", x$seed) else ""))
  fmt <- function(m, v) {
    pct_metrics <- c("sensitivity", "specificity", "ppv_study_controls",
                     "npv", "ppv_population")
    if (m %in% pct_metrics) sprintf("%.2f%%", 100 * v) else format(round(v), big.mark = ",")
  }
  for (i in seq_len(nrow(x$ci))) {
    m <- x$ci$metric[i]
    cat(sprintf("  %-19s %s  (%s to %s)\n", m,
                fmt(m, x$ci$estimate[i]), fmt(m, x$ci$lower[i]), fmt(m, x$ci$upper[i])))
  }
  has_pop <- "ppv_population" %in% x$ci$metric
  cat("\nPPV is the absolute risk of the reaction in a carrier (test-positive).\n")
  cat("  ppv_study_controls: anchored on this study's own controls (their specificity).\n")
  if (has_pop) {
    cat("  ppv_population:     anchored on the population allele carriage (generalizable).\n")
  }
  if (length(x$warnings)) {
    cat("\nNotes:\n")
    for (w in x$warnings) cat(" -", w, "\n")
  }
  invisible(x)
}

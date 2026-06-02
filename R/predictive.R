#' Predictive values, number needed to test, and carrier risk
#'
#' Translates an HLA risk-allele association into clinical-utility metrics for
#' use as a screening test. The positive predictive value equals the absolute
#' risk of the reaction in an exposed carrier.
#'
#' Identities used (carrier/dominant coding; rare reaction, so OR approx RR):
#' \deqn{PPV = sensitivity \times incidence / carriage}
#' \deqn{NPV = 1 - (1 - sensitivity) \times incidence / (1 - carriage)}
#' \deqn{NNT = 1 / (sensitivity \times incidence)}
#' The PPV ceiling (at 100\% sensitivity) is `incidence / carriage`, and the
#' carrier risk is enriched over the population incidence by the factor
#' `sensitivity / carriage`.
#'
#' @param sensitivity Proportion of cases carrying the allele (test sensitivity),
#'   in \[0, 1\].
#' @param carriage Population carriage frequency of the allele among the exposed
#'   population, in (0, 1).
#' @param incidence Incidence of the reaction among exposed individuals, as a
#'   proportion (e.g. `1/1000`). See [one_in()].
#' @param per_flagged Denominator for the natural-frequency breakdown
#'   (default 1000 flagged carriers).
#'
#' @return An object of class `hla_translation`: a list with `ppv`, `npv`,
#'   `nnt`, `carrier_risk`, `ppv_ceiling`, `enrichment_vs_pop`, `miss_fraction`,
#'   `noncarrier_risk`, `react_per_flagged`, `unnecessary_per_flagged`, and
#'   `inputs`. A `print()` method gives a readable summary.
#'
#' @examples
#' # Lamotrigine, HLA-A*32:01
#' ppv_npv_nnt(sensitivity = 0.414, carriage = 0.06, incidence = one_in(1000))
#' @export
ppv_npv_nnt <- function(sensitivity, carriage, incidence, per_flagged = 1000) {
  stopifnot(
    is.numeric(sensitivity), sensitivity >= 0, sensitivity <= 1,
    is.numeric(carriage),    carriage    >  0, carriage    <  1,
    is.numeric(incidence),   incidence   >= 0, incidence   <  1,
    is.numeric(per_flagged), per_flagged >  0
  )

  ppv         <- min(1, sensitivity * incidence / carriage)
  npv         <- 1 - (1 - sensitivity) * incidence / (1 - carriage)
  nnt         <- if (sensitivity * incidence > 0) 1 / (sensitivity * incidence) else Inf
  ppv_ceiling <- min(1, incidence / carriage)
  enrichment  <- sensitivity / carriage
  react       <- ppv * per_flagged
  resid_nc    <- (1 - sensitivity) * incidence / (1 - carriage)

  structure(list(
    ppv                     = ppv,
    npv                     = npv,
    nnt                     = nnt,
    carrier_risk            = ppv,
    ppv_ceiling             = ppv_ceiling,
    enrichment_vs_pop       = enrichment,
    miss_fraction           = 1 - sensitivity,
    noncarrier_risk         = resid_nc,
    react_per_flagged       = react,
    unnecessary_per_flagged = per_flagged - react,
    inputs = list(sensitivity = sensitivity, carriage = carriage,
                  incidence = incidence, per_flagged = per_flagged)
  ), class = "hla_translation")
}

#' Convert a "1 in N" frequency to a proportion
#'
#' @param N The denominator of a "1 in N" rate.
#' @return The proportion `1 / N`.
#' @examples
#' one_in(1000)
#' @export
one_in <- function(N) {
  stopifnot(is.numeric(N), N > 0)
  1 / N
}

#' Positive predictive value from an odds ratio
#'
#' Cross-check on [ppv_npv_nnt()] that derives PPV from the odds ratio rather
#' than the observed case carrier fraction, assuming a rare reaction
#' (OR approx RR).
#'
#' @inheritParams ppv_npv_nnt
#' @param or Odds ratio for carriers versus non-carriers.
#' @return The positive predictive value (numeric scalar).
#' @examples
#' ppv_from_or(or = 16.4, carriage = 0.06, incidence = one_in(1000))
#' @export
ppv_from_or <- function(or, carriage, incidence) {
  stopifnot(or > 0, carriage > 0, carriage < 1, incidence >= 0, incidence < 1)
  incidence * or / (carriage * or + (1 - carriage))
}

#' Positive predictive value from sensitivity, specificity, and incidence
#'
#' Prevalence-adjusted PPV using the test's specificity as the false-positive
#' anchor: \deqn{PPV = (Se \times I) / (Se \times I + (1 - Sp)(1 - I)).}
#' This is the estimator used by Manson et al. (2020,
#' doi:10.3389/fphar.2020.567048) and the appropriate form when specificity is
#' measured against tolerant controls rather than assumed equal to
#' `1 - carriage`. It coincides with [ppv_npv_nnt()]'s carriage-based PPV only
#' when `specificity == 1 - carriage`; when tolerant-control specificity departs
#' from population carriage (common in matched HLA designs) the two differ, and
#' this specificity-based form is the one to report against a tolerant-control
#' study.
#'
#' @param sensitivity Test sensitivity, in \[0, 1\].
#' @param specificity Test specificity (e.g. from tolerant controls), in \[0, 1\].
#' @param incidence Incidence of the reaction among exposed, as a proportion.
#' @return The positive predictive value (numeric scalar).
#' @examples
#' # Chung 2004 carbamazepine B*15:02: Se=1.0, Sp=0.970, incidence 0.0025
#' ppv_specificity(sensitivity = 1.0, specificity = 0.970, incidence = 0.0025)
#' @export
ppv_specificity <- function(sensitivity, specificity, incidence) {
  stopifnot(sensitivity >= 0, sensitivity <= 1,
            specificity >= 0, specificity <= 1,
            incidence   >= 0, incidence   <= 1)
  (sensitivity * incidence) /
    (sensitivity * incidence + (1 - specificity) * (1 - incidence))
}

#' @export
print.hla_translation <- function(x, ...) {
  ip <- x$inputs
  cat(sprintf("Inputs : sensitivity = %.1f%%,  carriage = %.2f%%,  incidence = %s (1 in %s)\n",
              100 * ip$sensitivity, 100 * ip$carriage,
              format(ip$incidence, scientific = FALSE),
              format(round(1 / ip$incidence), big.mark = ",")))
  cat(sprintf("PPV    : %.2f%%   (= risk in an exposed carrier; %.1fx the population rate)\n",
              100 * x$ppv, x$enrichment_vs_pop))
  cat(sprintf("NPV    : %.3f%%  (1 in %s non-carriers still react)\n",
              100 * x$npv,
              if (x$noncarrier_risk > 0) format(round(1 / x$noncarrier_risk), big.mark = ",") else "Inf"))
  cat(sprintf("NNT    : %s   (test this many, withhold from carriers, to prevent one case)\n",
              format(round(x$nnt), big.mark = ",")))
  cat(sprintf("Of %s flagged carriers: %s react, %s avoid the trigger unnecessarily\n",
              format(ip$per_flagged, big.mark = ","),
              format(round(x$react_per_flagged), big.mark = ","),
              format(round(x$unnecessary_per_flagged), big.mark = ",")))
  cat(sprintf("PPV ceiling (100%% sensitivity) = %.2f%% ; test misses %.0f%% of cases\n",
              100 * x$ppv_ceiling, 100 * x$miss_fraction))
  invisible(x)
}

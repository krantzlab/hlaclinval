#' Predictive values with confidence intervals from a 2x2 table
#'
#' Computes sensitivity, specificity, PPV, NPV and number needed to test with
#' confidence intervals from raw case-control carrier counts, by delegating to
#' [epiR::epi.tests()]. Confidence intervals require cell counts (not summary
#' rates), so this is the count-based companion to [ppv_npv_nnt()], which takes
#' population rates and returns point estimates.
#'
#' The 2x2 is laid out as test (carrier status) by outcome (case/control):
#' \tabular{lcc}{
#'             \tab case      \tab control     \cr
#'   carrier   \tab `cases_carrier`   \tab `controls_carrier`   \cr
#'   non-carrier \tab `cases_noncarrier` \tab `controls_noncarrier` \cr
#' }
#'
#' Note: in a case-control sample, PPV/NPV from `epi.tests` reflect the sampled
#' case:control ratio, not the population. To obtain population PPV/NPV at a
#' realistic incidence, supply `prevalence` (incidence among exposed); the
#' function then also returns prevalence-adjusted PPV/NPV using the
#' sensitivity/specificity point estimates.
#'
#' @param cases_carrier,cases_noncarrier Carrier and non-carrier counts in cases.
#' @param controls_carrier,controls_noncarrier Carrier and non-carrier counts in
#'   controls.
#' @param conf.level Confidence level (default 0.95).
#' @param prevalence Optional incidence of the reaction among exposed; if given,
#'   prevalence-adjusted PPV/NPV are added.
#'
#' @return A list with `epir` (the full [epiR::epi.tests()] object, including
#'   sensitivity, specificity, PPV, NPV and their CIs) and, when `prevalence` is
#'   supplied, `ppv_adj` and `npv_adj`.
#'
#' @examples
#' \dontrun{
#'   # Konvinse vancomycin DRESS: 19/23 cases, 0/46 controls carried A*32:01
#'   predictive_ci(cases_carrier = 19, cases_noncarrier = 4,
#'                 controls_carrier = 0, controls_noncarrier = 46)
#' }
#' @export
predictive_ci <- function(cases_carrier, cases_noncarrier,
                          controls_carrier, controls_noncarrier,
                          conf.level = 0.95, prevalence = NULL) {
  if (!requireNamespace("epiR", quietly = TRUE)) {
    stop("Package 'epiR' is required for confidence intervals. ",
         "Install it with install.packages('epiR'), or use ppv_npv_nnt() ",
         "for point estimates.", call. = FALSE)
  }
  # epi.tests expects: test +/- (rows) by outcome +/- (cols), as
  # c(TP, FP, FN, TN) filled by row.
  dat <- as.table(matrix(
    c(cases_carrier,    controls_carrier,
      cases_noncarrier, controls_noncarrier),
    nrow = 2, byrow = TRUE,
    dimnames = list(Test = c("carrier", "noncarrier"),
                    Outcome = c("case", "control"))
  ))
  res <- epiR::epi.tests(dat, conf.level = conf.level)

  out <- list(epir = res)

  if (!is.null(prevalence)) {
    se <- res$detail$est[res$detail$statistic == "se"]
    sp <- res$detail$est[res$detail$statistic == "sp"]
    out$ppv_adj <- (se * prevalence) /
      (se * prevalence + (1 - sp) * (1 - prevalence))
    out$npv_adj <- (sp * (1 - prevalence)) /
      ((1 - se) * prevalence + sp * (1 - prevalence))
  }
  out
}

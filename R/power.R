#' Exact power for a case-control carrier-status association
#'
#' Computes the power of the two-sided Fisher exact test to detect a
#' carrier-status association, by enumerating the exact rejection region over
#' all possible 2x2 tables and summing the joint probability under two
#' independent binomials (carriers among cases and among controls). This is
#' exact (no simulation, no normal approximation) and conservative at small
#' sample sizes.
#'
#' @param n_case Number of cases.
#' @param carriage Control / population carriage frequency (exposure prevalence),
#'   in (0, 1).
#' @param or Odds ratio to detect.
#' @param alpha Significance threshold. For a Bonferroni correction over `m`
#'   tests, pass `0.05 / m`.
#' @param ratio Controls per case (ignored if `n_ctrl` is supplied).
#' @param n_ctrl Number of controls. Overrides `ratio` when supplied.
#' @param reject Optional precomputed rejection matrix from
#'   [build_reject()]; supply it to avoid recomputation across a grid of ORs.
#'
#' @return Power (numeric scalar in \[0, 1\]).
#' @examples
#' # 29 cases, 2:1 controls, 4.1% carriage, single test
#' fisher_power(n_case = 29, carriage = 0.041, or = 16.4, ratio = 2)
#' @export
fisher_power <- function(n_case, carriage, or, alpha = 0.05,
                         ratio = 1, n_ctrl = NULL, reject = NULL) {
  stopifnot(n_case >= 1, carriage > 0, carriage < 1, or > 0, alpha > 0, alpha < 1)
  if (is.null(n_ctrl)) n_ctrl <- round(n_case * ratio)
  if (is.null(reject)) reject <- build_reject(n_case, n_ctrl, alpha)
  p1 <- case_prop(or, carriage)
  w1 <- stats::dbinom(0:n_case, n_case, p1)
  w2 <- stats::dbinom(0:n_ctrl, n_ctrl, carriage)
  sum(outer(w1, w2) * reject)
}

#' Minimum detectable odds ratio at a target power
#'
#' Finds the smallest odds ratio detectable at `power` using exact Fisher power
#' ([fisher_power()]). Sample size is treated as fixed (as in retrospective
#' registry/biobank studies); this is the recommended quantity for
#' STROBE/STREGA item 10 rather than post-hoc observed power.
#'
#' @inheritParams fisher_power
#' @param power Target power (default 0.80).
#' @param or_max Upper search bound for the odds ratio (default 80).
#' @return The minimum detectable odds ratio, or `NA` if `power` is not
#'   reachable below `or_max`.
#' @examples
#' # Minimum detectable OR: 29 cases, 2:1, 4.1% carriage, alpha 0.05
#' min_detectable_or(n_case = 29, carriage = 0.041, ratio = 2)
#'
#' # With Bonferroni over 100 alleles
#' min_detectable_or(n_case = 29, carriage = 0.041, ratio = 2, alpha = 0.05 / 100)
#' @export
min_detectable_or <- function(n_case, carriage, alpha = 0.05,
                              ratio = 1, n_ctrl = NULL,
                              power = 0.80, or_max = 80) {
  stopifnot(power > 0, power < 1, or_max > 1)
  if (is.null(n_ctrl)) n_ctrl <- round(n_case * ratio)
  reject <- build_reject(n_case, n_ctrl, alpha)
  f <- function(or) {
    fisher_power(n_case, carriage, or, alpha, n_ctrl = n_ctrl, reject = reject) - power
  }
  if (f(or_max)  < 0) return(NA_real_)
  if (f(1.0001) >= 0) return(1)
  stats::uniroot(f, c(1.0001, or_max), tol = 1e-3)$root
}

#' Expected case carriage from an odds ratio
#'
#' @inheritParams fisher_power
#' @return The expected carrier frequency among cases.
#' @export
case_prop <- function(or, carriage) {
  (or * carriage) / (1 - carriage + or * carriage)
}

#' Exact rejection region for the two-sided Fisher exact test
#'
#' Builds the logical matrix `R[x1 + 1, x2 + 1]` indicating whether a 2x2 table
#' with `x1` carriers among cases and `x2` among controls is rejected at
#' `alpha`. Depends only on the sample sizes and `alpha` (not the odds ratio),
#' so it can be built once and reused across an OR grid.
#'
#' @param n_case Number of cases.
#' @param n_ctrl Number of controls.
#' @param alpha Significance threshold.
#' @return A logical matrix of dimension `(n_case + 1) x (n_ctrl + 1)`.
#' @export
build_reject <- function(n_case, n_ctrl, alpha) {
  R <- matrix(FALSE, n_case + 1L, n_ctrl + 1L)
  for (x1 in 0:n_case) {
    for (x2 in 0:n_ctrl) {
      R[x1 + 1L, x2 + 1L] <- fisher_two_sided_p(x1, x2, n_case, n_ctrl) <= alpha
    }
  }
  R
}

#' Two-sided Fisher exact p-value for one 2x2 table
#'
#' Computed via the hypergeometric distribution conditional on the margins.
#'
#' @param x1 Carriers among cases.
#' @param x2 Carriers among controls.
#' @param n_case Number of cases.
#' @param n_ctrl Number of controls.
#' @return The two-sided Fisher exact p-value.
#' @export
fisher_two_sided_p <- function(x1, x2, n_case, n_ctrl) {
  N  <- n_case + n_ctrl
  m1 <- x1 + x2
  lo <- max(0, n_case - (N - m1)); hi <- min(n_case, m1)
  ks  <- lo:hi
  pr  <- stats::dhyper(ks, m1, N - m1, n_case)
  obs <- stats::dhyper(x1, m1, N - m1, n_case)
  sum(pr[pr <= obs * (1 + 1e-7)])
}

#' Minor allele frequency from carriage (dominant model, HWE)
#'
#' Convenience for cross-checking against packages parameterised by minor allele
#' frequency (e.g. genpwr). Under Hardy-Weinberg equilibrium, carriage equals
#' `1 - (1 - MAF)^2`, so `MAF = 1 - sqrt(1 - carriage)`.
#'
#' @param carriage Population carriage frequency, in (0, 1).
#' @return The implied minor allele frequency.
#' @examples
#' maf_from_carriage(0.06)
#' @export
maf_from_carriage <- function(carriage) {
  stopifnot(carriage > 0, carriage < 1)
  1 - sqrt(1 - carriage)
}

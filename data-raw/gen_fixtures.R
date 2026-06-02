#!/usr/bin/env Rscript
###############################################################################
## gen_fixtures.R
## ---------------------------------------------------------------------------
## Single source of truth for the web calculators. Runs the `hlaclinval` package
## over a grid of inputs and writes a JSON fixtures file. The browser widgets
## (PPV/NPV/NNT and exact-Fisher power) test their JavaScript implementations
## against these values, so the JS and the R package cannot silently diverge.
##
## Usage:
##   Rscript gen_fixtures.R [output_path]
## Default output: hlaclinval_fixtures.json
##
## No external packages required (hand-rolled JSON writer, base-R only).
###############################################################################

suppressMessages(library(hlaclinval))

out_path <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(out_path)) out_path[1] else "hlaclinval_fixtures.json"

## ---- tiny base-R JSON writer (numbers, strings, named lists, arrays) -------
json_num <- function(x) {
  if (!is.finite(x)) return("null")
  formatC(x, format = "g", digits = 12)
}
to_json <- function(x, indent = 0) {
  pad  <- strrep("  ", indent)
  padc <- strrep("  ", indent + 1)
  if (is.null(x)) return("null")
  if (is.character(x) && length(x) == 1) return(paste0("\"", x, "\""))
  if (is.numeric(x) && length(x) == 1)  return(json_num(x))
  if (is.logical(x) && length(x) == 1)  return(if (x) "true" else "false")
  if (is.list(x) && !is.null(names(x))) {
    items <- vapply(names(x), function(nm)
      paste0(padc, "\"", nm, "\": ", to_json(x[[nm]], indent + 1)), character(1))
    return(paste0("{\n", paste(items, collapse = ",\n"), "\n", pad, "}"))
  }
  # unnamed list or atomic vector -> array
  xs <- if (is.list(x)) x else as.list(x)
  items <- vapply(xs, function(e) paste0(padc, to_json(e, indent + 1)), character(1))
  paste0("[\n", paste(items, collapse = ",\n"), "\n", pad, "]")
}

## ---- grids -----------------------------------------------------------------
round6 <- function(x) if (is.finite(x)) signif(x, 12) else x

predictive_grid <- function() {
  sens <- c(0.10, 0.414, 0.50, 0.826, 0.98, 1.00)
  carr <- c(0.02, 0.06, 0.063, 0.09, 0.168)
  Ns   <- c(37, 50, 68, 417, 1000, 10000, 1e5)
  cases <- list()
  for (s in sens) for (f in carr) for (N in Ns) {
    r <- ppv_npv_nnt(s, f, one_in(N))
    cases[[length(cases) + 1]] <- list(
      sensitivity = s, carriage = f, incidence_one_in = N,
      ppv          = round6(r$ppv),
      npv          = round6(r$npv),
      nnt          = round6(r$nnt),
      carrier_risk = round6(r$carrier_risk),
      ppv_ceiling  = round6(r$ppv_ceiling)
    )
  }
  cases
}

power_grid <- function() {
  combos <- list(
    list(n_case = 29,  carriage = 0.041, ratio = 2, alpha = 0.05),
    list(n_case = 29,  carriage = 0.041, ratio = 2, alpha = 0.05/100),
    list(n_case = 23,  carriage = 0.063, ratio = 2, alpha = 0.05),
    list(n_case = 50,  carriage = 0.10,  ratio = 1, alpha = 0.05),
    list(n_case = 100, carriage = 0.06,  ratio = 2, alpha = 0.05/20000)
  )
  ors <- c(2, 5, 10, 16.4, 30)
  cases <- list()
  for (cc in combos) {
    n_ctrl <- round(cc$n_case * cc$ratio)
    reject <- build_reject(cc$n_case, n_ctrl, cc$alpha)
    pw <- lapply(ors, function(o) list(
      or    = o,
      power = round6(fisher_power(cc$n_case, cc$carriage, o, cc$alpha,
                                  n_ctrl = n_ctrl, reject = reject))
    ))
    md <- min_detectable_or(cc$n_case, cc$carriage, cc$alpha,
                            n_ctrl = n_ctrl, power = 0.80)
    cases[[length(cases) + 1]] <- list(
      n_case = cc$n_case, n_ctrl = n_ctrl, carriage = cc$carriage,
      alpha = cc$alpha, mdor_80 = round6(md), power_at_or = pw
    )
  }
  cases
}

fixtures <- list(
  package      = "hlaclinval",
  generated    = as.character(Sys.Date()),
  r_version    = as.character(getRversion()),
  description  = "Expected values for the web calculators; regenerate with gen_fixtures.R",
  predictive   = predictive_grid(),
  power        = power_grid()
)

writeLines(to_json(fixtures), out_path)
cat(sprintf("Wrote %d predictive + %d power fixtures to %s\n",
            length(fixtures$predictive), length(fixtures$power), out_path))

## ---- CI fixtures (R-side regression; not used by the JS widgets) -----------
## Appended block: regenerate writes a companion file with seeded CI bounds.
if (requireNamespace("hlaclinval", quietly = TRUE)) {
  ci <- suppressWarnings(hlaclinval::ppv_npv_nnt_ci(
    cases_carrier = 19, cases_n = 23, controls_carrier = 0, controls_n = 46,
    incidence = c(1/100, 1/40), population_carriers = 3418, population_n = 54249,
    seed = 1, n_boot = 20000))
  ci_fix <- list(
    scenario = "konvinse_vanco_seed1_nboot20000",
    seed = 1, n_boot = 20000,
    nnt_lower = round(ci$ci$lower[ci$ci$metric == "nnt"], 3),
    nnt_upper = round(ci$ci$upper[ci$ci$metric == "nnt"], 3),
    ppv_population_lower = round(ci$ci$lower[ci$ci$metric == "ppv_population"], 6),
    ppv_population_upper = round(ci$ci$upper[ci$ci$metric == "ppv_population"], 6)
  )
  writeLines(to_json(ci_fix), sub("\\.json$", "_ci.json", out_path))
  cat("Wrote CI fixture to", sub("\\.json$", "_ci.json", out_path), "\n")
}

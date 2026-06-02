# hlaclinval

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20514792.svg)](https://doi.org/10.5281/zenodo.20514792)

Clinical-translation metrics for HLA association studies. One package covering both the predictive-value layer and the study-design layer, because a SCAR analysis uses them in sequence: *was I powered to detect this association, and if I found it, what is its clinical utility?*

## Install

```r
# install.packages("remotes")
remotes::install_github("krantzlab/hlaclinval")
```

Core functions are base-R only. `testthat` is used for the test suite (Suggests).

## Predictive value, NPV, and NNT

The positive predictive value of an HLA risk allele used as a screening test is
the absolute risk of the reaction in an exposed carrier:

```r
library(hlaclinval)

# Lamotrigine, HLA-A*32:01
ppv_npv_nnt(sensitivity = 0.414, carriage = 0.06, incidence = one_in(1000))
#> PPV    : 0.69%   (= risk in an exposed carrier; 6.9x the population rate)
#> NPV    : 99.938% (1 in 1,604 non-carriers still react)
#> NNT    : 2,415
#> PPV ceiling (100% sensitivity) = 1.67% ; test misses 59% of cases
```

- `sensitivity` — proportion of cases carrying the allele
- `carriage` — population carriage frequency among the exposed
- `incidence` — incidence of the reaction among exposed (`one_in(1000)` = 0.001)

`ppv_from_or()` gives an odds-ratio-based cross-check. Note that PPV is bounded
above by `incidence / carriage`, so a low PPV for a rare reaction is structural,
not a defect of the marker.

## Exact power and minimum detectable OR

Power is computed exactly from the Fisher exact test (no simulation, no normal
approximation), appropriate for the small, sparse samples typical of SCAR
genetics.

```r
# 29 cases, 2:1 controls, 4.1% carriage
min_detectable_or(n_case = 29, carriage = 0.041, ratio = 2)       #> ~8.2
fisher_power(n_case = 29, carriage = 0.041, or = 16.4, ratio = 2) #> ~0.99

# Bonferroni over 100 alleles
min_detectable_or(n_case = 29, carriage = 0.041, ratio = 2, alpha = 0.05 / 100) #> ~17.5
```

For retrospective registry/biobank studies, report sample size as fixed by case
availability and supplement with the minimum detectable OR (STROBE/STREGA item
10), rather than post-hoc observed power.

To cross-check against packages parameterised by allele frequency (e.g.
`genpwr`), convert with `maf_from_carriage()` (MAF = 1 - sqrt(1 - carriage),
dominant model under HWE).

## Confidence intervals

Point estimates from small SCAR samples can mislead. `ppv_npv_nnt_ci()` returns
intervals that propagate the sampling uncertainty in sensitivity and specificity
*and* the (usually dominant) uncertainty in incidence:

```r
# Konvinse vancomycin: 19/23 cases, 0/46 controls, incidence uncertain
ppv_npv_nnt_ci(cases_carrier = 19, cases_n = 23,
               controls_carrier = 0, controls_n = 46,
               incidence = c(1/100, 1/40),
               population_carriers = 3418, population_n = 54249, seed = 1)
#> nnt           77  (48 to 130)
#> ppv_population  20.73%  (12.24% to 33.28%)   # brackets the measured 19.2% carrier risk
```

Sensitivity and specificity are drawn from Jeffreys Beta distributions (a
smoothed bootstrap that stays well-behaved at zero cells, unlike resampling
counts); Se/Sp themselves are reported with exact Clopper-Pearson intervals;
NNT uses a percentile interval because it is a skewed ratio. Incidence is an
explicit uncertain input -- pass `incidence = c(lo, hi)` (log-uniform) or a
`function(n)` sampler. A fixed scalar incidence triggers a warning that the
resulting intervals understate total uncertainty.

## Manuscript outputs

The package returns objects, never files, and depends on no plotting or table
library -- it produces the numbers and leaves formatting to you.

```r
# Tidy data frame -> pipe into flextable / gt / knitr::kable for a Word table
as.data.frame(ci_result)

# Power-justification table across designs
power_grid(n_case = c(20, 29, 50), carriage = 0.041, ratio = 2,
           alpha = c(0.05, 0.05/100))

# Signature figures (base-R graphics)
plot_carrier_risk_vs_incidence(sensitivity = 0.414, carriage = 0.06,
                               mark_incidence = 1/1000)   # utility is set by incidence
plot_power_curve(n_case = 29, carriage = 0.041, ratio = 2, mark_or = 16.4)
```

## Validation

The test suite reproduces published HLA-drug associations:

| Example | sensitivity | carriage | incidence | NNT | source |
|---|---|---|---|---|---|
| Carbamazepine B\*15:02 | 0.98 | 16.8% | 1 in 417 | ~426 | CPIC ~461 |
| Vancomycin A\*32:01 | 0.826 | 6.3% | 1 in 68* | ~82 | Konvinse 19.2% carrier risk |
| Lamotrigine A\*32:01 | 0.414 | 6.0% | 1 in 1,000 | ~2,415 | this study |

\* The vancomycin incidence is back-derived from Konvinse's directly measured
19.2% carrier risk, shown only to reproduce that figure; carbamazepine B\*15:02
(with an independently published incidence) is the clean forward validation.

## Caveats

Carrier (dominant) coding and a rare reaction (OR approx RR) are assumed.
Incidence is typically the least precise input and PPV/NNT scale with it, so
report results across a plausible incidence range rather than a single point.

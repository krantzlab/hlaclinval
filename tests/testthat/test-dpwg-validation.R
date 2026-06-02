## Validation against Manson, Swen & Guchelaar (2020), "Diagnostic Test Criteria
## for HLA Genotyping to Prevent Drug Hypersensitivity Reactions", Front.
## Pharmacol. 11:567048 (doi:10.3389/fphar.2020.567048).
##
## The review computed, per primary study,
##   NNG = 1 / (sensitivity * incidence)
##   PPV = (Se * I) / (Se * I + (1 - Sp) * (1 - I))   [tolerant-control specificity]
## We reproduce both against published Table values:
##   - NNG validates the sensitivity/incidence machinery shared with min_detectable_or()
##   - specificity-based PPV validates ppv_specificity()
## Values transcribed verbatim from Tables 2-13.

dpwg <- data.frame(
  label = c("ABC Mallal2008 immuno (T2)", "ABC Saag white immuno (T2)",
            "ALLO Chiu2012 SCAR (T3)", "ALLO Cheng2015 SCAR (T3)",
            "ALLO Cao2012 SJS/TEN (T5)", "FLUCLOX Daly2009 DILI (T6)",
            "CBZ Chung2004 SJS/TEN (T7)", "CBZ Cheung2013 SJS/TEN (T7)",
            "CBZ Hsiao2014 SJS/TEN (T7)", "LTG Shi2017 SJS/TEN (T13)",
            "LTG Koomdee2017 SJS/TEN (T13)"),
  se  = c(1.000, 1.000, 1.000, 0.946, 1.000, 0.843, 1.000, 0.923, 0.884, 0.227, 0.250),
  sp  = c(0.969, 0.960, 0.867, 0.880, 0.889, 0.938, 0.970, 0.881, 0.928, 0.814, 0.880),
  inc = c(0.027, 0.026, 0.0021, 0.0021, 0.0016, 0.000085, 0.0025, 0.0025, 0.0025, 0.001, 0.001),
  ppv = c(0.479, 0.406, 0.0155, 0.0163, 0.0142, 0.0011, 0.0781, 0.0191, 0.0297, 0.0012, 0.0021),
  nng = c(37, 38, 476, 504, 625, 13953, 400, 433, 453, 4400, 4000),
  stringsAsFactors = FALSE
)

test_that("DPWG review: NNG reproduced (1 / [sensitivity * incidence])", {
  for (i in seq_len(nrow(dpwg))) {
    r <- ppv_npv_nnt(dpwg$se[i], carriage = 0.10, incidence = dpwg$inc[i])
    # published NNG are rounded (often to nearest 100); allow 1% or 5, whichever larger
    expect_equal(round(r$nnt), dpwg$nng[i], tolerance = 0.02, info = dpwg$label[i])
  }
})

test_that("DPWG review: specificity-based PPV reproduced", {
  for (i in seq_len(nrow(dpwg))) {
    p <- ppv_specificity(dpwg$se[i], dpwg$sp[i], dpwg$inc[i])
    expect_equal(p, dpwg$ppv[i], tolerance = 0.05,
                 info = dpwg$label[i])
  }
})

test_that("carriage-based PPV is the rare-disease limit of the specificity form", {
  # Se*I/f equals (Se*I)/(Se*I+(1-Sp)(1-I)) with Sp=1-f only as I -> 0.
  # Agreement should tighten as incidence shrinks.
  d <- function(I) abs(ppv_npv_nnt(0.8, 0.06, I)$ppv - ppv_specificity(0.8, 1 - 0.06, I))
  expect_lt(d(1e-3), 5e-4)     # SCAR-scale: negligible
  expect_lt(d(1e-5), 1e-7)     # ultra-rare: converges
  expect_gt(d(1e-2), 1e-3)     # common phenotype: forms genuinely differ
})

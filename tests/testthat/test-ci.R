## Tests for the Monte Carlo CI family (ppv_npv_nnt_ci).

test_that("point estimates are deterministic and match closed forms", {
  r <- suppressWarnings(ppv_npv_nnt_ci(19, 23, 0, 46, incidence = 1/68, seed = 1))
  expect_equal(unname(r$estimate["sensitivity"]), 19/23, tolerance = 1e-9)
  expect_equal(unname(r$estimate["specificity"]), 1.0,  tolerance = 1e-9)
  # specificity-based PPV with Sp=1 pins at 1
  expect_equal(unname(r$estimate["ppv_study_controls"]), 1.0, tolerance = 1e-9)
  expect_equal(unname(r$estimate["nnt"]), 1 / (19/23 * 1/68), tolerance = 1e-6)
})

test_that("intervals are ordered and bracket the estimate", {
  r <- suppressWarnings(ppv_npv_nnt_ci(24, 26, 16, 135, incidence = c(1/500, 1/250), seed = 42))
  for (i in seq_len(nrow(r$ci))) {
    expect_lte(r$ci$lower[i], r$ci$estimate[i] + 1e-9)
    expect_gte(r$ci$upper[i], r$ci$estimate[i] - 1e-9)
    expect_lte(r$ci$lower[i], r$ci$upper[i])
  }
})

test_that("seed makes results reproducible", {
  a <- suppressWarnings(ppv_npv_nnt_ci(24, 26, 16, 135, incidence = c(1/500,1/250), seed = 7))
  b <- suppressWarnings(ppv_npv_nnt_ci(24, 26, 16, 135, incidence = c(1/500,1/250), seed = 7))
  expect_identical(a$ci, b$ci)
})

test_that("fixed incidence warns and narrows PPV vs an incidence range", {
  expect_warning(
    withCallingHandlers(
      ppv_npv_nnt_ci(24, 26, 16, 135, incidence = 1/400, seed = 1),
      warning = function(w) if (grepl("Sparse data", conditionMessage(w)))
        invokeRestart("muffleWarning")
    ),
    "Incidence treated as fixed"
  )
  fixed <- suppressWarnings(ppv_npv_nnt_ci(24, 26, 16, 135, incidence = 1/400, seed = 1))
  rng   <- suppressWarnings(ppv_npv_nnt_ci(24, 26, 16, 135, incidence = c(1/800,1/200), seed = 1))
  w <- function(r) diff(as.numeric(r$ci[r$ci$metric == "ppv_study_controls", c("lower","upper")]))
  expect_lt(w(fixed), w(rng))
})

test_that("zero control cell warns and the carriage-based PPV stays realistic", {
  expect_warning(ppv_npv_nnt_ci(19, 23, 0, 46, incidence = c(1/100,1/40),
                                population_carriers = 3418, population_n = 54249, seed = 1),
                 "Sparse data")
  r <- suppressWarnings(ppv_npv_nnt_ci(19, 23, 0, 46, incidence = c(1/100,1/40),
                                       population_carriers = 3418, population_n = 54249,
                                       seed = 1, n_boot = 20000))
  # carriage-based PPV brackets Konvinse's measured 19.2% carrier risk
  cr <- r$ci[r$ci$metric == "ppv_population", ]
  expect_lt(cr$lower, 0.192); expect_gt(cr$upper, 0.192)
  # NNT interval brackets the published ~75 (seeded; allow MC slack)
  nn <- r$ci[r$ci$metric == "nnt", ]
  expect_lt(nn$lower, 75); expect_gt(nn$upper, 75)
  expect_equal(nn$lower, 47.8, tolerance = 4)
  expect_equal(nn$upper, 129,  tolerance = 8)
})

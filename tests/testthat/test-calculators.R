test_that("PPV/NPV/NNT match the lamotrigine A*32:01 example", {
  r <- ppv_npv_nnt(sensitivity = 0.414, carriage = 0.06, incidence = one_in(1000))
  expect_equal(r$ppv,          0.0069, tolerance = 1e-3)   # ~0.69%
  expect_equal(r$nnt,          2415,   tolerance = 1)       # ~2,415
  expect_equal(r$enrichment_vs_pop, 6.9, tolerance = 0.05)  # ~6.9x population rate
  expect_equal(r$ppv_ceiling,  1/0.06 * 1/1000, tolerance = 1e-6)
  expect_true(r$npv > 0.999)
})

test_that("Carbamazepine B*15:02 reproduces the published NNT (~461)", {
  # CPIC: ~461 patients need testing to prevent one CBZ-SJS/TEN case
  r <- ppv_npv_nnt(sensitivity = 0.98, carriage = 0.168, incidence = one_in(417))
  expect_equal(r$nnt, 461, tolerance = 40)   # exact 1/(0.98*1/417) = 425; CPIC ~461
  expect_true(r$ppv < 0.05)                  # low PPV, as published
  expect_true(r$npv > 0.999)
})

test_that("Vancomycin A*32:01 reproduces Konvinse's 19.2% carrier risk", {
  # Carrier risk (PPV) measured directly by Konvinse = 19.2%
  r <- ppv_npv_nnt(sensitivity = 0.826, carriage = 0.063, incidence = one_in(68))
  expect_equal(r$carrier_risk, 0.192, tolerance = 0.01)
  expect_equal(r$nnt, 82, tolerance = 5)      # ~82, vs published ~75
})

test_that("ppv_from_or agrees with the sensitivity-based PPV in range", {
  p_or  <- ppv_from_or(or = 16.4, carriage = 0.06, incidence = one_in(1000))
  expect_equal(p_or, 0.00852, tolerance = 1e-2)   # ~0.85%
})

test_that("one_in inverts a 1-in-N rate", {
  expect_equal(one_in(1000), 0.001)
  expect_equal(one_in(50), 0.02)
})

test_that("input validation rejects out-of-range arguments", {
  expect_error(ppv_npv_nnt(sensitivity = 1.5, carriage = 0.06, incidence = 0.001))
  expect_error(ppv_npv_nnt(sensitivity = 0.4, carriage = 0,    incidence = 0.001))
  expect_error(ppv_npv_nnt(sensitivity = 0.4, carriage = 0.06, incidence = 2))
})

test_that("fisher_power is monotone in OR and bounded in [0,1]", {
  p_lo <- fisher_power(n_case = 29, carriage = 0.041, or = 2,  ratio = 2)
  p_hi <- fisher_power(n_case = 29, carriage = 0.041, or = 20, ratio = 2)
  expect_gt(p_hi, p_lo)
  expect_gte(p_lo, 0); expect_lte(p_hi, 1)
})

test_that("observed lamotrigine OR is well powered; MDOR is sensible", {
  pw <- fisher_power(n_case = 29, carriage = 0.041, or = 16.4, ratio = 2, alpha = 0.05)
  expect_gt(pw, 0.90)                         # ~99%
  md <- min_detectable_or(n_case = 29, carriage = 0.041, ratio = 2, alpha = 0.05)
  expect_equal(md, 8.2, tolerance = 1.0)      # ~8.2 at 80% power
})

test_that("Bonferroni correction raises the minimum detectable OR", {
  md1   <- min_detectable_or(n_case = 29, carriage = 0.041, ratio = 2, alpha = 0.05)
  md100 <- min_detectable_or(n_case = 29, carriage = 0.041, ratio = 2, alpha = 0.05/100)
  expect_gt(md100, md1)
})

test_that("unreachable target power returns NA", {
  md <- min_detectable_or(n_case = 5, carriage = 0.04, ratio = 1,
                          alpha = 0.05/1e5, power = 0.80, or_max = 80)
  expect_true(is.na(md))
})

test_that("maf_from_carriage inverts the HWE carriage relation", {
  maf <- maf_from_carriage(0.06)
  expect_equal(1 - (1 - maf)^2, 0.06, tolerance = 1e-8)
})

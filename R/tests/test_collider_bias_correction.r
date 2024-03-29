library(testthat)
library(vroom)
source("R/util.r")
source("R/collider_bias.r")

test_that("collider_bias.correct_for_collider_bias works well", {
  collider_bias_results_file <- tempfile(fileext = ".tsv")
  slopehunter_file <- tempfile(fileext = ".tsv.gz")
  harmonised_results <- tempfile(fileext = ".tsv.gz")

  suppressWarnings(
      conduct_collider_bias_analysis("R/tests/data/test_data_small.tsv.gz",
                                   "R/tests/data/test_data_small.tsv.gz",
                                   "R/tests/data/clumped_snps.tsv.gz",
                                   collider_bias_results_file,
                                   harmonised_results,
                                   slopehunter_file,
                                   p_value_thresholds = c(0.001)
    )
  )

  expect_true(file.exists(collider_bias_results_file))
  result <- vroom::vroom(collider_bias_results_file)
  expect_equal(nrow(result), 3)
  expect_true(all(result$BETA >= 1))

  expect_true(file.exists(harmonised_results))
  expect_true(file.exists(slopehunter_file))

  slopehunter_result <- vroom::vroom(slopehunter_file)
  expect_true(all(c("P", "BETA", "SE") %in% colnames(slopehunter_result)))
})

library(testthat)
library(vroom)
source("R/util.r")

test_that("ldsc.sh ensure that parsing the file works well", {
  collider_bias_results_file <- tempfile(fileext = ".tsv")
  slopehunter_file <- tempfile(fileext = ".tsv.gz")
  harmonised_results <- tempfile(fileext = ".tsv.gz")


  expect_true(file.exists(collider_bias_results_file))
  result <- vroom::vroom(collider_bias_results_file)
  expect_equal(nrow(result), 3)
  expect_true(all(result$BETA >= 1))

  expect_true(file.exists(harmonised_results))
  expect_true(file.exists(slopehunter_file))

  slopehunter_result <- vroom::vroom(slopehunter_file)
  expect_true(all(c("P", "BETA", "SE") %in% colnames(slopehunter_result)))
})

library(testthat)
library(vroom)
source("r_scripts/functions/util.r")
source("r_scripts/functions/collider_bias.r")

test_that("collider_bias.correct_for_collider_bias works well", {
  collider_bias_results_file <- "/tmp/collider_bias_results.tsv"
  slopehunter_file <- "/tmp/slopehunter.tsv.gz"
  harmonised_results <- "/tmp/harmonised_results.tsv.gz"
  correct_for_collider_bias("r_scripts/tests/data/test_data_small.tsv.gz",
                            "r_scripts/tests/data/test_data_small.tsv.gz",
                            "r_scripts/tests/data/clumped_snps.tsv.gz",
                            collider_bias_results_file,
                            harmonised_results,
                            slopehunter_file,
                            p_value_thresholds = c(0.001)
  )

  expect_true(file.exists(collider_bias_results_file))
  result <- vroom::vroom(collider_bias_results_file)
  expect_equal(nrow(result), 3)

  expect_true(file.exists(slopehunter_file))
  expect_true(file.exists(harmonised_results))
})

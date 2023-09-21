library(testthat)
source("R/mr.r")
source("R/util.r")
source("R/graphs.r")

pqtl_top_hits_dir <- "R/tests/data"
mr_results <- tempfile(fileext = ".tsv")

test_that("mr.perform_mr_on_pqtl_datasets runs cis only mr against a gwas", {
  perform_mr_on_pqtl_datasets("R/tests/data/test_data_small.tsv.gz", mr_results)

  expect_true(file.exists(mr_results))
  results <- vroom::vroom(mr_results)
  expect_equal(nrow(results), 12)
})

test_that("mr.compare_interesting_mr_results creates a grouped forest plot", {
  forest_plot_output_file <- tempfile(fileext = ".png")
  compare_interesting_mr_results(mr_results, forest_plot_output_file)

  expect_true(file.exists(forest_plot_output_file))
  expect_gt(file.info(forest_plot_output_file)$size, 10000)
})

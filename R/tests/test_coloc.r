library(testthat)
library(vroom)
source("R/coloc.r")
source("R/util.r")
source("R/gwas_formatting.r")

test_that("coloc.coloc_analysis returns a successfull coloc analysis", {
  output_file <- tempfile(fileext = ".tsv")

  result <- coloc_analysis("R/tests/data/test_data_small.tsv.gz",
                           "R/tests/data/test_data_small.tsv.gz",
                           output_file,
                           exposure_name = "exposure",
                           chr = 1,
                           bp = 1232132,
                           range = 5000000
  )
  expect_equal(nrow(result), 1)
  expect_true(all(result$h0 > 0.95 ))
  expect_true(all(result$h1 < 0.01 ))
})

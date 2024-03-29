library(testthat)
source("R/coloc.r")
source("R/util.r")
source("R/gwas_formatting.r")

test_that("coloc.coloc_analysis returns a successfull coloc analysis", {
  output_file <- tempfile(fileext = ".tsv")
  gwas <- vroom::vroom("R/tests/data/test_data_small.tsv.gz")

  #suppressing warning around coloc not using great data from the GWAS
  result <- suppressWarnings(coloc_analysis(gwas,
                           gwas,
                           exposure_name = "exposure",
                           chr = 1,
                           bp = 1232132,
                           range = 5000000
  ))
  expect_equal(nrow(result), 1)
  expect_true(all(result$h0 > 0.95 ))
  expect_true(all(result$h1 < 0.01 ))
})

test_that("coloc.run_coloc_on_list_of_datasets returns a file out coloc results", {
  output_file <- tempfile(fileext = ".tsv")
  gwas <- "R/tests/data/test_data_small.tsv.gz"

  #suppressing warning around coloc not using great data from the GWAS
  result <- suppressWarnings(run_coloc_on_list_of_datasets(c(gwas),
                                            c(gwas),
                                            exposure_name = c("exposure"),
                                            chr = c(1),
                                            bp = c(1232132),
                                            range = c(5000000)
  ))
  expect_equal(nrow(result), 1)
  expect_true(all(result$h0 > 0.95 ))
  expect_true(all(result$h1 < 0.01 ))
})

library(testthat)
library(vroom)
source("r_scripts/functions/coloc.r")
source("r_scripts/functions/util.r")
source("r_scripts/functions/gwas_formatting.r")

test_that("coloc.coloc_analysis returns a successfull coloc analysis", {
  output_file <- "/tmp/coloc_result.tsv"
  file.remove(output_file)

  result <- coloc_analysis("r_scripts/tests/data/test_data_small.tsv.gz",
                 "r_scripts/tests/data/test_data_small.tsv.gz",
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

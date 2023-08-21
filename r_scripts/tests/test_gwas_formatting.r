library(testthat)
library(vroom)
source("r_scripts/functions/gwas_formatting.r")
source("r_scripts/functions/util.r")

test_gwas <- vroom::vroom("r_scripts/tests/data/test_data_small.tsv.gz")

test_that("gwas_formatting.standardise_gwas standardises a gwas", {
  output_file <- "/tmp/test_data_standardised.tsv.gz"
  standardise_gwas("r_scripts/tests/data/test_data_small.tsv.gz", output_file, "default")
  result <- vroom::vroom(output_file)

  expect_equal(nrow(result), nrow(test_gwas))
  expect_true(all(result$EA < result$OA))
  expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
})

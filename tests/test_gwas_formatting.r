library(testthat)
library(vroom)

test_gwas <- vroom::vroom("data/test_data_small.tsv.gz")
gwas_num_rows <- nrow(test_gwas)

that_that("gwas_formatting.standardise_gwas standardises a gwas", {
  source("../r_scripts/functions/gwas_formatting.r")
  output_file <- "/tmp/test_data_standardised.tsv.gz"
  standardise_gwas("data/test_data_small.tsv.gz", output_file, "default")
  result <- vroom:vroom(output_file)

  expect_equal(nrow(result), gwas_num_rows)
  expect_true(all(result$EA < result$OA))
  expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
  expect_true(is.null(result$RSID))
})

# that_that("gwas_formatting.standardise_gwas standardises a gwas and populates RSID", {
#   source("../r_scripts/functions/gwas_formatting.r")
#   output_file <- "/tmp/test_data_standardised.tsv.gz"
#   standardise_gwas("data/test_data_small.tsv.gz", output_file, "default", populate_rsid = T)
#   result <- vroom:vroom(output_file)
#
#   expect_equal(nrow(result), gwas_num_rows)
#   expect_true(all(result$EA < result$OA))
#   expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
#   expect_true(is.null(result$RSID))
# })

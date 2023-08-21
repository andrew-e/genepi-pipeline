library(testthat)
library(vroom)
source("r_scripts/functions/util.r")
source("r_scripts/functions/gwas_formatting.r")
source("r_scripts/functions/gwas_comparisons.r")

test_that("gwas_comparison.compare_two_gwases_from_clumped_hits works well", {
  result <- compare_two_gwases_from_clumped_hits("r_scripts/tests/data/test_data_small.tsv.gz",
                                       "r_scripts/tests/data/test_data_small.tsv.gz",
                                       "r_scripts/tests/data/clumped_snps.tsv.gz"
  )

  expect_equal(nrow(result$res), 4)
  expect_true(all(result$res$VALUE > 28))
  expect_equal(nrow(result$variants), 29)
})

test_that("gwas_comparison.run_cross_section_of_gwas_comparisons works well", {
  input_gwases <- c("r_scripts/tests/data/test_data_small.tsv.gz", "r_scripts/tests/data/test_data_small.tsv.gz", "r_scripts/tests/data/test_data_small.tsv.gz")
  input_clumps <- c("r_scripts/tests/data/clumped_snps.tsv.gz", "r_scripts/tests/data/clumped_snps.tsv.gz", "r_scripts/tests/data/clumped_snps.tsv.gz")
  output_file <- "/tmp/expected_vs_observed_results.tsv"
  variants_output_file <- "/tmp/expected_vs_observed_variants.tsv"

  compare_replication_across_all_gwas_permutations(input_gwases, input_clumps, output_file, variants_output_file)

  expect_true(file.exists(output_file))
  results <- vroom::vroom(output_file)
  expect_equal(nrow(results), 4 * length(input_gwases))
  expect_true(all(results$VALUE > 28))

  expect_true(file.exists(variants_output_file))
  variants <- vroom::vroom(variants_output_file)
  expect_equal(nrow(variants), 29 * length(input_gwases))
})

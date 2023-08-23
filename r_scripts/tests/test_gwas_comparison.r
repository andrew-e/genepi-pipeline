library(testthat)
library(vroom)
source("r_scripts/functions/util.r")
source("r_scripts/functions/graphs.r")
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

test_that("gwas_comparison.compare_heterogeneity_across_ancestries works well", {
  input_gwases <- c("r_scripts/tests/data/test_data_small.tsv.gz", "r_scripts/tests/data/test_data_small.tsv.gz", "r_scripts/tests/data/test_data_small.tsv.gz")
  input_clumps <- c("r_scripts/tests/data/clumped_snps.tsv.gz", "r_scripts/tests/data/clumped_snps.tsv.gz", "r_scripts/tests/data/clumped_snps.tsv.gz")
  input_ancestries <- c("EUR", "AFR", "SAS")

  heterogeneity_file <- "/tmp/heterogeneity_scores.tsv"
  first_plot <- "/tmp/heterogeneity_first_plot.png"
  second_plot <- "/tmp/heterogeneity_second_plot.png"

  compare_heterogeneity_across_ancestries(input_gwases, input_clumps, input_ancestries, heterogeneity_file, first_plot, second_plot)

  expect_true(file.exists(heterogeneity_file))
  results <- vroom::vroom(heterogeneity_file)
  print(nrow(results))
  #expect_equal(nrow(results), 4 * length(input_gwases))

  expect_true(file.exists(first_plot))
  expect_true(file.exists(second_plot))
})

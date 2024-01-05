library(testthat)
source("R/gwas_formatting.r")
source("R/data_conversions.r")
source("R/util.r")

test_that("gwas_formatting.standardise_gwas standardises a gwas", {
  test_gwas <- vroom::vroom("R/tests/data/test_data_small.tsv.gz")
  output_file  <- tempfile(fileext = ".tsv.gz")

  standardise_gwas("R/tests/data/test_data_small.tsv.gz", output_file, "default")
  result <- vroom::vroom(output_file)

  expect_equal(nrow(result), nrow(test_gwas))
  expect_true(all(result$EA < result$OA))
  expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
})

test_that("gwas_formatting.standardise_gwas with bespoke_column_map standardises a gwas", {
  output_file  <- tempfile(fileext = ".tsv.gz")
  map <- "SNP=MARKER,CHR=CHR,BP=BP,EA=A0,OA=A1,EAF=A0FREQ,P=P,BETA=BETA,SE=SE,OR=OR,OR_LB=OR_LB,OR_UB=OR_UB,RSID=RSID"
  bespoke_column_map <- parse_gwas_input_column_maps(map)
  bespoke_column_map <- split_string_into_named_list(bespoke_column_map)

  standardise_gwas("R/tests/data/test_data_tiny.tsv.gz", output_file, bespoke_column_map = bespoke_column_map)
  result <- vroom::vroom(output_file)

  expect_equal(nrow(result), 12)
  expect_true(all(result$EA < result$OA))
  expect_true(all(grep("\\d+:\\d+_\\w+_\\w+", result$SNP)))
})

test_that("gwas_formatting.convert_beta_to_or and back returns the same results", {
  original_gwas <- vroom::vroom("R/tests/data/test_data_small.tsv.gz")

  gwas <- convert_beta_to_or(original_gwas)
  gwas <- convert_or_to_beta(gwas)

  floating_point_tolerance <- 1e-10
  expect_true(all(abs(gwas$BETA - original_gwas$BETA) < floating_point_tolerance))
  expect_true(all(abs(gwas$SE- original_gwas$SE) < floating_point_tolerance))
})

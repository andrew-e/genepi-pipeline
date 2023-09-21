library(testthat)
library(vroom)
source("R/gwas_formatting.r")
source("R/util.r")

test_that("gwas_formatting.standardise_gwas standardises a gwas", {
  gwas <- "R/tests/data/test_data_small.tsv.gz"
  snp <- "19:12436574_A_G"
  result <- vroom_snps(gwas, snp, "10:100392738_C_T")
  print(result)
  expect_true(result[result$SNP == snp, ]$SNP == snp)
})

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

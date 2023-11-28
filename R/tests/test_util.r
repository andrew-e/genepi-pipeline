library(testthat)
library(vroom)
source("R/util.r")

test_that("util.get_file_or_dataframe returns a data frame if given a file", {
  file <- "R/tests/data/test_data_small.tsv.gz"
  snps <- c("19:12436574_A_G", "10:100392738_C_T")
  result <- get_file_or_dataframe(file, snps=snps)
  expect_true(all(result$SNP %in% snps))
  expect_equal(nrow(result), 2)
})

test_that("util.get_file_or_dataframe returns a data frame if given a file", {
  file <- "R/tests/data/test_data_small.tsv.gz"
  snps <- c("19:12436574_A_G", "10:100392738_C_T")
  columns <- c("SNP")
  result <- get_file_or_dataframe(file, columns=columns, snps=snps)
  expect_true(all(result$SNP %in% snps))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 1)
})

test_that("util.get_file_or_dataframe returns a data frame if given a file", {
  df <- vroom::vroom("R/tests/data/test_data_small.tsv.gz")
  snps <- c("19:12436574_A_G", "10:100392738_C_T")
  result <- get_file_or_dataframe(df, snps=snps)
  expect_true(all(result$SNP %in% snps))
  expect_equal(nrow(result), 2)
})

test_that("util.get_file_or_dataframe returns a data frame if given a file", {
  df <- vroom::vroom("R/tests/data/test_data_small.tsv.gz")
  snps <- c("19:12436574_A_G", "10:100392738_C_T")
  columns <- c("SNP")
  result <- get_file_or_dataframe(df, columns=columns, snps=snps)
  expect_true(all(result$SNP %in% snps))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 1)
})


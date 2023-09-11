library(testthat)
source("r_scripts/functions/util.r")

test_that( "r_markdown test rmd works", {
  input_string <- "test=/tmp/mr_results.tsv miami_plot=/tmp/mr_forest.png"
  rmd_file <- "/tmp/test.html"

  input_params <- split_string_into_named_list(input_string)
  create_html_from_rmd("r_scripts/markdown/test.rmd", input_params, rmd_file)

  expect_true(file.exists(rmd_file))
})

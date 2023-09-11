library(testthat)
source("r_scripts/functions/graphs.r")
source("r_scripts/functions/util.r")

test_that("coloc.coloc_analysis returns a successfull coloc analysis", {
  miami_plot_file <- "/tmp/miami_plot.png"
  file.remove(miami_plot_file)

  result <- miami_plot("r_scripts/tests/data/test_data_small.tsv.gz",
                       "r_scripts/tests/data/test_data_small.tsv.gz",
                       miami_plot_file,
                       title = "ooyooo",
                       chr = 1,
                       bp = 5232132,
                       range = 1000000
  )
  expect_true(file.exists(miami_plot_file))
  expect_gt(file.info(miami_plot_file)$size, 10000)
})
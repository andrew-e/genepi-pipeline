library(testthat)
source("R/graphs.r")
source("R/util.r")

test_that("graphs.miami_plot given a specific range, creates a png file ", {
  miami_plot_file <- tempfile(fileext = ".png")

  result <- miami_plot("R/tests/data/test_data_small.tsv.gz",
                       "R/tests/data/test_data_small.tsv.gz",
                       miami_plot_file,
                       title = "ooyooo",
                       chr = 1,
                       bp = 5232132,
                       range = 1000000
  )
  expect_true(file.exists(miami_plot_file))
  expect_gt(file.info(miami_plot_file)$size, 10000)
})
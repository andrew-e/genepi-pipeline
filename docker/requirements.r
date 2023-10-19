install.packages(c("R.utils", "argparser", "coloc", "devtools", "MendelianRandomization",
                   "eeptools", "gmp", "qqman", "corrplot", "broom", "rmarkdown", "ggrepel",
                   "conflicted", "nloptr", "Cairo", "plotrix", "forestplot", "dplyr", "BiocManager"),
                 repos = "http://cran.us.r-project.org")

install.packages('https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_23.8.tar.gz', repos = NULL, type = 'source')
.libPaths( c( .libPaths(), "/usr/lib/R/site-library") )


devtools::install_github(c("Osmahmoud/SlopeHunter",
                           "phenoscanner/phenoscanner",
                           "MRCIEU/ieugwasr",
                           "MRCIEU/TwoSampleMR",
                           "suchestoncampbelllab/gwasurvivr",
                           "yixuan/prettydoc")
)

BiocManager::install("GENESIS")

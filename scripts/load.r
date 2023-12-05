#maybe look at changing this to creating a package: https://stackoverflow.com/a/48094346/6104011
source("../R/util.r")
source("../R/constants.r")
lapply(list.files("../R/", full.names = T, pattern="\\.r$"), source)
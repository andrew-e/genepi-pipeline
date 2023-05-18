library(parallel)

#mcapply to do typical 'lapply()' things put all in parallel
detectCores()
r <- mclapply(1:10, function(i) {
    Sys.sleep(10)  ## Do nothing for 10 seconds
}, mc.cores = 10)  ## Split this job across 10 cores

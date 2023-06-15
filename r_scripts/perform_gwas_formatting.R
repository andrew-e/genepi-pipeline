library(readr)
library(vroom)
library(dplyr)
library(TwoSampleMR)

source("functions/gwas_formatting.r")

output_location <- ".../add/path/" # TODO - from argparse


data_lookup_file <- paste0("meta/gwas_summary_statistics.csv")
data_lookup<-read_csv(data_lookup_file)


for (i in 1:nrow(data_lookup)){
  
  trait <- data_lookup[i,]
  
  custom_cols <- trait %>% select(trait, starts_with("custom"))
  
  out <- 
    read_and_format(file_gwas = trait$original_file,
                  data_version= trait$source,
                  custom_file_columns=custom_cols)
  
  vroom_write(out, paste0(output_location, trait$trait, "_tidy_GWAS.txt.gz"))
  
  
}



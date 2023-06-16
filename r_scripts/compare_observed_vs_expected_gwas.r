source("home/r_scripts/load.r")
library(argparser, quietly=TRUE)

parser <- arg_parser("Comparison of expected vs observed replication rates")
parser <- add_argument(parser, "--gwas_filenames", help="Comma separated list of GWAS filenames", type="character")
parser <- add_argument(parser, "--clumped_filenames", help="Comma separated list of clumped SNP filenames", type="character")
parser <- add_argument(parser, "--result", help="Filename of results", type="character")

args <- parse_args(parser)

gwas_filenames <- split_string_into_list(args.gwas_filenames)
clumped_filenames<- split_string_into_list(args.cluimped_filenames)

for (i in seq_len(length(gwas_filenames))) {
    for (j in i:length(gwas_filenames)) {
        j = j+1
        if (j <= length(gwas_filenames)) {
            expected_vs_observed_replication(gwas_filenames[i], gwas_filenames[j], clumped_filenames[i], clumped_filenames[j])
        }
    }
}

#forest_plot <- function(sslist, snp)
#{
#    tibble(
#        beta = sapply(sslist, \(x) x$bhat[snp]),
#        se = sapply(sslist, \(x) x$se[snp]),
#        label = names(sslist)
#    ) %>%
#    ggplot(., aes(x=beta, y=label)) +
#    geom_point() +
#    geom_errorbarh(aes(xmin=beta-se*1.96, xmax=beta+se*1.96), height=0) +
#    geom_vline(xintercept=0, linetype="dotted") +
#    labs(x="beta", y="population")
#}
#
#alpha = 0.05/nrow(gwas_1$snp)
#
#list <- c("first_gwas", "second_gwas", "third_gwas", "fourth_gwas")
#
#index <- which(first_gwas$pval < 5e-8)
#o_first_second <- expected_vs_observed_replication(first_gwas$bhat[index], second_gwas$bhat[index], third_gwas$se[index], fourth_gwas$se[index], 0.05)
#
#png("expected_observed_fail_first_second_forest_plot.png", width = 4, height = 4, units = 'in', res = 300)
#forest_plot(list(eur=first_gwas, afr=second_gwas, eas=third_gwas, sas=fourth_gwas), index[which(o_first_second$variants$sign_fail)[3]])
#dev.off()

#TODO: a sort of specific grouped forest plot (maybe comparing ancestries?)
# forest_plot <- function(sslist, snp)
# {
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
# }
#
# alpha = 0.05/nrow(gwas_1$snp)
#
# list <- c("first_gwas", "second_gwas", "third_gwas", "fourth_gwas")
#
# index <- which(first_gwas$pval < 5e-8)
# o_first_second <- expected_vs_observed_replication(first_gwas$bhat[index], second_gwas$bhat[index], third_gwas$se[index], fourth_gwas$se[index], 0.05)
#
# png("expected_observed_fail_first_second_forest_plot.png", width = 4, height = 4, units = 'in', res = 300)
# forest_plot(list(eur=first_gwas, afr=second_gwas, eas=third_gwas, sas=fourth_gwas), index[which(o_first_second$variants$sign_fail)[3]])
# dev.off()

forest_plot <- function(table, title, output_file) {
  library(ggplot2, quietly = T)
  create_dir_for_files(output_file)

  if (!all(c("BETA", "SE") %in% names(table))) {
    stop("data frame needs to have BETA and SE columns")
  }

  first_column_name <- colnames(table)[1]
  table$BETA <- as.numeric(table$BETA)
  table$SE <- as.numeric(table$SE)

  table$LL <- table$BETA - (1.96 * table$SE)
  table$UL <- table$BETA + (1.96 * table$SE)

  ggplot(table,
         aes(y = .data[[first_column_name]], x = BETA, xmin = LL, xmax = UL, color = .data[[first_column_name]]) ) +
    geom_pointrange(cex = 1) +
    geom_vline(xintercept = 0) +
    theme(legend.position = "none") +
    ggtitle(title)

  ggsave(output_file)
}

grouped_forest_plot <- function(table, title, group_column, output_file, p_value_column = NA, q_stat_column = NA) {
  library(ggplot2, quietly = T)
  if (!("BETA" %in% names(table)) || !("SE" %in% names(table))) {
    stop("data frame needs to have BETA and SE named columns")
  }

  first_column_name <- colnames(table)[1]

  table$BETA <- as.numeric(table$BETA)
  table$SE <- as.numeric(table$SE)

  table$LL <- table$BETA - (1.96 * table$SE)
  table$UL <- table$BETA + (1.96 * table$SE)

  plot_thing <- ggplot(table,
         aes(y = if(!is.na(q_stat_column)) paste0(.data[[first_column_name]], "\n Q-stat=", .data[[q_stat_column]]) else .data[[first_column_name]],
             x = BETA,
             xmin = LL,
             xmax = UL,
             color=.data[[group_column]],
             fill=.data[[group_column]])) +
    ylab(first_column_name) +
    geom_vline(xintercept = 0) +
    geom_pointrange(cex = 1, fatten = 2.5, position=position_dodge(width = 0.5)) +
    theme(legend.position = "bottom") +
    ggtitle(title)

  if (!is.na(p_value_column)) {
    plot_thing + geom_text(aes(label = paste("P=", pval_adjusted), group = .data[[group_column]]), position = position_dodge(width = 0.5))
  }
  else if (!is.na(q_stat_column)) {
    plot_thing + geom_text(aes(x = 1.5, label = .data[[q_stat_column]]))
  }


  ggsave(output_file)
}

#' manhattan_and_qq: produce manhattan and qq plot from a GWAS file
#'
#' @param gwas_filename: a file of a gwas that includes CHR, CP, P, and SNP
#' @param name: name of plots to be saved (and named as a header in graph)
#' @param save_dir: defaults to 'scratch/results'
#' @return 2 plots: one manhattan plot and one QQ plot (with lambda included)
manhattan_and_qq <- function(gwas_filename, manhattan_filename, qq_filename, include_qq = T) {
  library(qqman, quietly = T)
  library(vroom, quietly = T)
  create_dir_for_files(manhattan_filename, qq_filename)

  manhattan_columns <- c("SNP", "CHR", "BP", "P")
  gwas <- vroom::vroom(gwas_filename, col_select = manhattan_columns)
  gwas <- gwas[complete.cases(gwas), ]

  png(manhattan_filename, width = 1500, height = 500)
  qqman::manhattan(gwas, main = "Manhattan plot of GWAS")
  dev.off()

  if (include_qq) {
    png(qq_filename, width = 500, height = 500)

    qqman::qq(gwas$P, main = "Q-Q plot of GWAS p-values")
    lambda <- median(qchisq(1 - gwas$P, 1)) / qchisq(0.5, 1)
    text(0.5, 4, paste("lambda", "=", signif(lambda, digits = 3)))

    dev.off()
  }
}


#' miami_plot: produce miami plot of GWAS data from two gwases
#'
#' @param gwas_dataframe: a dataframe that includes CHR, CP, P, and SNP
#' @param name: name of plots to be saved (and named as a header in graph)
#' @param save_dir: defaults to 'scratch/results'
#' 
miami_plot <- function(first_gwas_filename,
                       second_gwas_filename,
                       miami_plot_file,
                       title = "Comparing GWASes",
                       chr = NA,
                       bp = NA,
                       range = NA) {
  library(qqman, quietly=T)
  library(vroom, quietly=T)
  create_dir_for_files(miami_plot_file)

  show_specific_region <- !is.na(chr) & !is.na(bp) & !is.na(range)

  manhattan_columns <- c("SNP", "CHR", "BP", "P")
  first_gwas <- vroom::vroom(first_gwas_filename, col_select = manhattan_columns)
  second_gwas <- vroom::vroom(second_gwas_filename, col_select = manhattan_columns)
  first_gwas <- first_gwas[complete.cases(first_gwas), ]
  second_gwas <- second_gwas[complete.cases(second_gwas), ]

  png_width <- 1500
  png_height <- 800

  top_ylim <- max(-log10(second_gwas$P))
  x_range <- NULL
  x_lab <- ""
  if (show_specific_region) {
    png_width <- 900
    range <- floor(range / 2)
    x_range <- c(bp - range, bp + range)
    x_lab <- paste("Chromosome", chr)

    first_gwas <- subset(first_gwas, CHR == chr)
    second_gwas <- gwas_region(second_gwas, chr, bp, range)
    top_ylim <-  max(-log10(second_gwas$P))
  }

  png(miami_plot_file, width = png_width, height = png_height)
  par(mfrow = c(2, 1))
  par(mar = c(0, 5, 3, 3))

  if (show_specific_region) {
    qqman::manhattan(first_gwas, main = title, xlim = x_range)
  }
  else {
    qqman::manhattan(first_gwas, main = title)
  }

  par(mar = c(5, 5, 3, 3))
  qqman::manhattan(second_gwas, ylim = c(top_ylim, 0), xlab = x_lab, xaxt = "n")

  #if (show_specific_region) {
  #  text(0.5, 4, paste("# of SNPS Compared =", nrow(second_gwas)))
  #}
  dev.off()
}

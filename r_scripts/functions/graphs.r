forest_plot <- function(table, title, output_file) {
  library(ggplot2, quietly = T)

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
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(output_file, limitsize = F)
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
    scale_colour_brewer(type="qual") +
    geom_vline(xintercept = 0) +
    geom_pointrange(cex = 1, fatten = 2, position=position_dodge(width = 0.5)) +
    theme(legend.position = "bottom") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))

  if (!is.na(p_value_column)) {
    plot_thing + geom_text(aes(label = paste("P=", signif(.data[[p_value_column]]), digits=3), group = .data[[group_column]]), position = position_dodge(width = 0.5))
  }
  else if (!is.na(q_stat_column)) {
    plot_thing + geom_text(aes(x = 1.5, label = .data[[q_stat_column]]))
  }

  ggsave(output_file, width = 2000, units = "px")
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

    first_gwas <- gwas_region(first_gwas, chr, bp, range)
    second_gwas <- gwas_region(second_gwas, chr, bp, range)
    top_ylim <-  max(-log10(second_gwas$P))
  }
  else {
    first_chrs <- sort(unique(first_gwas$CHR))
    second_chrs <- sort(unique(second_gwas$CHR))
    shared_chrs <- sort(intersect(first_chrs, second_chrs))

    if (!(all(shared_chrs == first_chrs) & all(shared_chrs == second_gwas))) {
      first_gwas <- subset(first_gwas, CHR %in% shared_chrs)
      second_gwas <- subset(second_gwas, CHR %in% shared_chrs)
    }
  }
  print(nrow(first_gwas))
  print(nrow(second_gwas))

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
  if (show_specific_region) {
    qqman::manhattan(second_gwas, ylim = c(top_ylim, 0), xlim = x_range, xlab = x_lab, xaxt = "n")
  }
  else {
    qqman::manhattan(second_gwas, ylim = c(top_ylim, 0), xlab = x_lab, xaxt = "n")
  }

  dev.off()
}

plot_heritability_contribution_per_ancestry <- function(heterogeneity_results_qj, output_file) {
  library(ggplot2)
  library(tidyr)
  graph_width <- max(1000, nrow(heterogeneity_results_qj) * 100)
  plot <- tidyr::gather(as.data.frame(heterogeneity_results_qj), "key", "value", -SNP)

  ggplot(plot, aes(x=SNP, y=-log10(value))) +
    geom_point(aes(colour=key)) +
    scale_colour_brewer(type="qual") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Contribution to heterogeneity score broken down by population") +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(output_file, width = graph_width, units = "px")
}

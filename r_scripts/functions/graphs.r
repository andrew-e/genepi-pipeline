source("/home/r_scripts/functions/util.r")

forest_plot <- function(table,
                        name,
                        save_location = "scratch/results/") {

    if (!("BETA" %in% names(table)) || !("SE" %in% names(table))) {
        stop("data frame needs to have BETA and SE named columns")
    }

    first_column_name <- colnames(table)[1]
    file_to_save <- paste0(save_location, name, ".png")

    table$BETA <- as.numeric(table$BETA)
    table$SE <- as.numeric(table$SE)

    table$LL <- table$BETA - (1.96 * table$SE)
    table$UL <- table$BETA + (1.96 * table$SE)

    ggplot(table,
           aes(y = .data[[first_column_name]], x = BETA, xmin = LL, xmax = UL, color = .data[[first_column_name]])) +
           geom_pointrange(cex=1) +
           geom_vline(xintercept=0) +
           theme(legend.position="none")
    ggsave(file_to_save)
}

#' manhattan_and_qq: produce manhattan and qq plot from a GWAS file
#'
#' @param gwas_filename: a file of a gwas that includes CHR, CP, P, and SNP
#' @param name: name of plots to be saved (and named as a header in graph) 
#' @param save_dir: defaults to 'scratch/results' 
#' @return 2 plots: one manhattan plot and one QQ plot (with lambda included)
manhattan_and_qq <- function(gwas_filename, manhattan_filename, qq_filename, include_qq=T, columns=list()) {
    library(qqman)
    create_dir_for_files(c(manhattan_filename, qq_filename))
    gwas <- data.table::fread(gwas_filename)

    chr <- if(!is.null(columns$CHR)) columns$CHR else "CHR"
    bp <- if(!is.null(columns$BP)) columns$BP else "BP"
    p <- if(!is.null(columns$P)) columns$P else "P"
    snp <- if(!is.null(columns$MARKER)) columns$MARKER else "MARKER"

    gwas[[chr]] <- as.numeric(gwas[[chr]])
    gwas[[bp]] <- as.numeric(gwas[[bp]])
    gwas[[p]] <- as.numeric(gwas[[p]])

    png(manhattan_filename, width=1500, height=500)
    qqman::manhattan(gwas,
              chr = chr,
              bp = bp,
              p = p,
              snp = snp,
              main = "Manhattan plot of GWAS")
    dev.off()

    if (include_qq) {
        png(qq_filename, width=500, height=500)

        qqman::qq(gwas[[p]], main = "Q-Q plot of GWAS p-values")

        lambda <- median(qchisq(1 - gwas[[p]], 1)) / qchisq(0.5,1)
        text(0.5, 4, paste("lambda", "=",  signif(lambda, digits = 3)))
        dev.off()
    }
}


#' miami_plot: produce miami plot of GWAS data
#'
#' @param gwas_dataframe: a dataframe that includes CHR, CP, P, and SNP
#' @param name: name of plots to be saved (and named as a header in graph)
#' @param save_dir: defaults to 'scratch/results'
#' @return 2 plots: one manhattan plot and one QQ plot (with lambda included)
#' @examples
miami_plot <- function(first_gwas_filename, second_gwas_filename, miami_filename, title="Comparing GWASes", columns=list()) {
    first_gwas <- data.table::fread(first_gwas_filename)
    second_gwas <- data.table::fread(second_gwas_filename)

    chr <- if(!is.null(columns$CHR)) columns$CHR else "CHR"
    bp <- if(!is.null(columns$BP)) columns$BP else "BP"
    p <- if(!is.null(columns$P)) columns$P else "P"
    snp <- if(!is.null(columns$MARKER)) columns$MARKER else "MARKER"

    first_gwas[[chr]] <- as.numeric(first_gwas[[chr]])
    first_gwas[[bp]] <- as.numeric(first_gwas[[bp]])
    first_gwas[[p]] <- as.numeric(first_gwas[[p]])
    first_gwas <- first_gwas[complete.cases(first_gwas),]

    second_gwas[[chr]] <- as.numeric(second_gwas[[chr]])
    second_gwas[[bp]] <- as.numeric(second_gwas[[bp]])
    second_gwas[[p]] <- as.numeric(second_gwas[[p]])
    second_gwas <- second_gwas[complete.cases(second_gwas),]

    png(miami_filename, width=1500, height=500)

    par(mfrow=c(2,1))
    par(mar=c(0,5,3,3))
    manhattan(first_gwas,
          ylim=c(0,10),
          chr = chr,
          bp = bp,
          p = p,
          snp = snp,
          main = title)

    par(mar=c(5,5,3,3))
    manhattan(second_gwas,
          ylim=c(10,0),
          chr = chr,
          bp = bp,
          p = p,
          snp = snp,
          xlab = "",
          xaxt = "n"
    )
    dev.off()
}

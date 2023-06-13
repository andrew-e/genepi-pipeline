source("/home/scripts/functions/util.r")

#' convert_odds_ratio_to_beta_estimate: Given an OR and lower and upper bounds,
#'   calculates the BETA, and SE.
#'   based on this answer: https://stats.stackexchange.com/a/327684
#'
#' @param df: dataframe with the following columns: OR, LB (lower bound), UB (upper bound)
#' @return df with new columns BETA and SE
#'
convert_odds_ratio_to_beta_estimate <- function(gwas_filename, columns=list()) {
    library(boot)

    gwas <- data.table::fread(gwas_filename)
    gwas$BETA <- log(gwas$OR)

    lower_bound <- boot::inv.logit(gwas$LB)
    upper_bound <- boot::inv.logit(gwas$UB)

    gwas$SE <- (gwas$UB - gwas$LB) / (2 * 1.95996)

    return(gwas)
}

#' expected_vs_observed_replication: compares two gwases by comparing significant/clumped list of SNPS between
#'   two different GWASes
#'
#' @param b_disc vector of clumped incidence hit effects
#' @param se_disc the standard errors for incidence effects
#' @param b_rep corresponding vector of associations in progression
#' @param se_rep standard errors of effects in progression
#' @param alpha p-value threshold to check for replication of incidence hits in progression (e.g. try 0.05 or 1e-5)
#'
expected_vs_observed_replication <- function(b_disc, b_rep, se_disc, se_rep, alpha = 0.05) {
    p_sign <- pnorm(-abs(b_disc) / se_disc) *
            pnorm(-abs(b_disc) / se_rep) +
            (
                (1 - pnorm(-abs(b_disc) / se_disc)) *
                (1 - pnorm(-abs(b_disc) / se_rep)
            )
    )

    p_sig <- pnorm(-abs(b_disc) / se_rep + qnorm(alpha / 2)) +
            (1 - pnorm(-abs(b_disc) / se_rep - qnorm(alpha / 2)))
    p_rep <- pnorm(abs(b_rep) / se_rep, lower.tail = FALSE)

    res <- tibble::tibble(
        nsnp = length(b_disc),
        metric = c("Sign", "Sign", "P-value", "P-value"),
        datum = c("Expected", "Observed", "Expected", "Observed"),
        value = c(
            sum(p_sign, na.rm=TRUE),
            sum(sign(b_disc) == sign(b_rep)),
            sum(p_sig, na.rm=TRUE),
            sum(p_rep < alpha, na.rm=TRUE)
        )
      )

    return(list(res = res, variants = dplyr::tibble(sig = p_sig, sign = p_sign)))
}

#' harmoise_gwases: takes a list of gwas filenames, and ensures that all effect allele, other allele,
#'  frequency, and beta values are aligned, then saves the files
#'
#' @param: gwas_filenames
harmonise_gwases <- function(gwas_filenames) {

}
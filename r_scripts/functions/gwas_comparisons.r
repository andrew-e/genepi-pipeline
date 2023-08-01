source("functions/util.r")

#' Expected vs observed replication rates
#'
#' @description For a set of effects that have discovery and replication betas and SEs, this function determines the extent to which the observed replication rate matches the expected replication rate.
#' The expected replication rate is based on the assumption that the replication dataset has the same effect sizes but that the power may be different (e.g. due to allele frequencies or sample sizes) and is reflected in the replication standard errors.
#' It assesses replication based on concordance of effect direction across discovery and replication, and p-values surpassing a user-specified p-value threshold.
#'
#' @param b_disc vector of clumped incidence hit effects
#' @param se_disc the standard errors for incidence effects
#' @param b_rep corresponding vector of associations in progression
#' @param se_rep standard errors of effects in progression
#' @param alpha p-value threshold to check for replication of incidence hits in progression (e.g. try 0.05 or 1e-5)

#TODO: Change this to input 2 gwas files and 2 clumped snp files?
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

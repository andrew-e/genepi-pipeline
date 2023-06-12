
#' correct_for_collider_bias: do a bunch of collider bias corrections.  Including:
#'  * SlopeHunter
#'  * Dudbridge Correction
#'  * IVW Correction
#'  * Otherss?
#'
#' @return 2 plots: one manhattan plot and one QQ plot (with lambda included)
#' @examples
correct_for_collider_bias <- function(incidence_gwas,
                                      subsequent_gwas,
                                      clumped_snps,
                                      save_as = F,
                                      p_value_threshold = 0.001,
                                      include_slopehunter = T,
                                      include_dudbridge = T,
                                      include_mr = T,
                                      produce_plots = F) {
    suppressPackageStartupMessages(library(SlopeHunter))

    incidence <- SlopeHunter::read_incidence(incidence_gwas,
                                snp_col = "MARKER",
                                effect_allele_col = "EA",
                                other_allele_col = "NONEA",
                                eaf_col = "EAF",
                                pval_col = "P",
                                beta_col = "BETA",
                                se_col = "SE",
                                chr_col = "CHR",
                                pos_col = "BP")

    incidence <- incidence[!(is.na(incidence$EA.incidence)|is.na(incidence$OA.incidence)),]
    subsequent_progression <- SlopeHunter::read_prognosis(subsequent_gwas,
                                             snp_col = "MARKER",
                                             effect_allele_col = "EA",
                                             other_allele_col = "NONEA",
                                             eaf_col = "EAF",
                                             pval_col = "P",
                                             beta_col = "BETA",
                                             se_col = "SE",
                                             chr_col = "CHR",
                                             pos_col = "BP")

    subsequent_progression <- subsequent_progression[!(is.na(subsequent_progression$EA.prognosis)|is.na(subsequent_progression$OA.prognosis)),]

    harmonised_effects <- SlopeHunter::harmonise_effects(incidence_dat = incidence,
                                            prognosis_dat = subsequent_progression,
                                            by.pos = FALSE,
                                            pos_cols = c("POS.incidence", "POS.prognosis"),
                                            snp_cols = c("SNP", "SNP"),
                                            beta_cols = c("BETA.incidence", "BETA.prognosis"),
                                            se_cols = c("SE.incidence", "SE.prognosis"),
                                            EA_cols = c("EA.incidence", "EA.prognosis"),
                                            OA_cols = c("OA.incidence", "OA.prognosis")
    )

    harmonised_effects <- harmonised_effects[!harmonised_effects$remove, ]
    harmonised_effects <- subset(harmonised_effects, remove==FALSE | (palindromic==TRUE & remove==TRUE))

    clumped_snps$MARKER <- tolower(clumped_snps$MARKER)
    pruned_harmonised_effects <- subset(harmonised_effects, (SNP %in% clumped_snps$MARKER) == TRUE)
    pruned_harmonised_effects <- subset(pruned_harmonised_effects, is.na(BETA.prognosis) == FALSE)
    pruned_harmonised_effects <- subset(pruned_harmonised_effects, BETA.prognosis < 10)

    comparing_results <- data.frame(METHOD = character(),
        BETA = numeric(),
        SE = numeric(),
        SNPS_USED = numeric()
    )

    if (include_slopehunter) {
        print("Starting SlopeHunter")

        column_name <- "SLOPEHUNTER"
        pruned_harmonised_effects_df <- data.frame(pruned_harmonised_effects)
        slopehunter_result <- SlopeHunter::hunt(dat = pruned_harmonised_effects_df,
                                   snp_col = "SNP",
                                   xbeta_col = "BETA.incidence",
                                   xse_col = "SE.incidence",
                                   xp_col = "PVAL.incidence",
                                   ybeta_col = "BETA.prognosis",
                                   yse_col = "SE.prognosis",
                                   yp_col = "PVAL.prognosis",
                                   xp_thresh = p_value_threshold
        )

        slopehunter_estimated_slope <- slopehunter_result$b
        slopehunter_estimated_slope_standard_error <- slopehunter_result$bse

        entropy <- slopehunter_result$entropy
        fit <- slopehunter_result$Fit
        hunted <- table(fit$clusters)[1]
        pleiotropic <- table(fit$clusters)[2]

        print("SlopeHunter slope and SE, and data")
        print(slopehunter_estimated_slope)
        print(slopehunter_estimated_slope_standard_error)
        print(paste("entropy: ", entropy))
        print(paste("hunted: ", hunted))
        print(paste("pleiotropic: ", pleiotropic))

        harmonised_effects <- adjust_gwas_data_from_weights(harmonised_effects, column_name, slopehunter_estimated_slope, slopehunter_estimated_slope_standard_error)

        if (produce_plots) {
            plot_name_prefix <- paste0(gsub(".*/(.+)\\.\\w*", "\\1", subsequent_gwas), "_joined_slopehunter")
            columns <- list(CHR="CHR.incidence", BP="POS.incidence", P=paste0("P_", column_name), MARKER="SNP")
            produce_plots(harmonised_effects, plot_name_prefix, columns=columns)

            comparing_results <- comparing_betas %>% add_row(METHOD = "SlopeHunter",
                                                           BETA = slopehunter_estimated_slope,
                                                           SE = slopehunter_estimated_slope_standard_error,
                                                           SNPS_USED = hunted)
        }
    }

    if (include_dudbridge) {
        suppressPackageStartupMessages(library(MendelianRandomization))
        print("Starting Dudbridge")

        column_name <- "DUDBRIDGE"

        weighted_dudbridge <- MendelianRandomization::mr_ivw(mr_input(bx = pruned_harmonised_effects$BETA.incidence,
                                              bxse = pruned_harmonised_effects$SE.incidence,
                                              by = pruned_harmonised_effects$BETA.prognosis,
                                              byse = pruned_harmonised_effects$SE.prognosis))

        pruned_harmonised_effects$weights <- 1 / pruned_harmonised_effects$SE.prognosis^2

        weighting <- (sum(pruned_harmonised_effects$weights * pruned_harmonised_effects$BETA.incidence^2)) /
          (
              (sum(pruned_harmonised_effects$weights * pruned_harmonised_effects$BETA.incidence ^ 2)) -
              (sum(pruned_harmonised_effects$weights * pruned_harmonised_effects$SE.incidence ^ 2))
          )

        dudbridge_estimated_slope <- weighted_dudbridge$Estimate * weighting
        dudbridge_estimated_slope_standard_error <- weighted_dudbridge$StdError * weighting

        print("Dudbridge slope and SE")
        print(dudbridge_estimated_slope)
        print(dudbridge_estimated_slope_standard_error)

        harmonised_effects <- adjust_gwas_data_from_weights(harmonised_effects, column_name, dudbridge_estimated_slope, dudbridge_estimated_slope_standard_error)

        if (produce_plots) {
            plot_name_prefix <- paste0(gsub(".*/(.+)\\.\\w*", "\\1", subsequent_gwas), "joined_dudbridge")

            columns <- list(CHR="CHR.incidence", BP="POS.incidence", P=paste0("P_", column_name), MARKER="SNP")
            produce_plots(harmonised_effects, plot_name_prefix, columns=columns)

            comparing_results <- comparing_betas %>% add_row(METHOD = "Dudbridge",
                                                           BETA = dudbridge_estimated_slope,
                                                           SE = dudbridge_estimated_slope_standard_error,
                                                           SNPS_USED = 0)

        }
    }


    if (produce_plots) {
        plot_name_prefix <- paste0(gsub(".*/(.+)\\.\\w*", "\\1", subsequent_gwas), "_joined_comparison")
        create_forest_plot(comparing_results, plot_name_prefix)
    }

    if (isTruthy(save_as)) {
        write.table(harmonised_effects, save_as, row.names=F, quote=F)
    } else {
        return(harmonised_effects)
    }
}

#' adjust_gwas_data_from_weights: Apply a slope correction weight (and SE) to an existing GWAS
#'   This can be used in conjuction with weights calculated to account for collider bias.
#'   Currently used to work on a GWAS result of Slopehunter.
#'
#' @param gwas: a dataframe that includes BETA.incidence, BETA.prognosis, SE.incidence, SE.prognosis
#' @param type_of_adjustment: string name of adjustment (eg. SLOPEHUNTER)
#' @param slope: number, slope of the correction
#' @param slope_standard_error: number, SE of the corrected slope
#' @return gwas with 3 additional columns BETA_{name}, SE_{name}, and P_{name}
bias.adjust_gwas_data_from_weights <- function(gwas, type_of_adjustment, slope, slope_standard_error) {
    adjusted_beta <- paste0("BETA_", type_of_adjustment)
    adjusted_se <- paste0("SE_", type_of_adjustment)
    adjusted_p <- paste0("P_", type_of_adjustment)

    gwas[[adjusted_beta]] <- ( gwas$BETA.prognosis - ( slope * gwas$BETA.incidence ) )

    gwas[[adjusted_se]] <- sqrt(
        (gwas$SE.prognosis * gwas$SE.prognosis) +
        ((slope * slope) * (gwas$SE.incidence * gwas$SE.incidence)) +
        ((gwas$BETA.incidence * gwas$BETA.incidence) * (slope_standard_error * slope_standard_error)) +
        ((gwas$SE.incidence * gwas$SE.incidence) * (slope_standard_error * slope_standard_error))
    )

    gwas[[adjusted_p]] <- pchisq(
        (gwas[[adjusted_beta]] / gwas[[paste0("SE_", type_of_adjustment)]])^2, 1, lower.tail = FALSE
    )

    return(gwas)
}

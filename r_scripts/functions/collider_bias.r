#' correct_for_collider_bias: do a bunch of collider bias corrections.  Including:
#'  * SlopeHunter
#'  * Dudbridge Correction
#'
#' @return 2 plots: one manhattan plot and one QQ plot (with lambda included)
#' @examples
correct_for_collider_bias <- function(incidence_gwas,
                                      subsequent_gwas,
                                      clumped_snps = c(),
                                      output_file,
                                      #TODO: maybe change p_value_theshold to list?
                                      p_value_threshold = 0.001,
                                      include_slopehunter = T,
                                      include_dudbridge = T) {
  suppressPackageStartupMessages(library(SlopeHunter))

  incidence <- SlopeHunter::read_incidence(incidence_gwas,
    snp_col = "SNP",
    effect_allele_col = "EA",
    other_allele_col = "OA",
    eaf_col = "EAF",
    pval_col = "P",
    beta_col = "BETA",
    se_col = "SE",
    chr_col = "CHR",
    pos_col = "BP"
  )

  incidence <- incidence[!(is.na(incidence$EA.incidence) | is.na(incidence$OA.incidence)), ]
  subsequent_progression <- SlopeHunter::read_prognosis(subsequent_gwas,
    snp_col = "SNP",
    effect_allele_col = "EA",
    other_allele_col = "NONEA",
    eaf_col = "EAF",
    pval_col = "P",
    beta_col = "BETA",
    se_col = "SE",
    chr_col = "CHR",
    pos_col = "BP"
  )

  subsequent_progression <- subsequent_progression[!(is.na(subsequent_progression$EA.prognosis) | is.na(subsequent_progression$OA.prognosis)), ]

  harmonised_effects <- SlopeHunter::harmonise_effects(
    incidence_dat = incidence,
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
  harmonised_effects <- subset(harmonised_effects, remove == FALSE | (palindromic == TRUE & remove == TRUE))

  pruned_harmonised_effects <- subset(harmonised_effects, (SNP %in% clumped_snps) == TRUE)
  pruned_harmonised_effects <- subset(pruned_harmonised_effects, is.na(BETA.prognosis) == FALSE)
  pruned_harmonised_effects <- subset(pruned_harmonised_effects, BETA.prognosis < 10)

  collider_bias_results <- data.frame(
    METHOD = character(),
    P_VAL_THRESHOLD = numeric(),
    BETA = numeric(),
    SE = numeric(),
    HUNTED = numeric(),
    PLEIOTROPIC = numeric(),
    ENTROPY = numeric()
  )

  if (include_slopehunter) {
    print("Starting SlopeHunter")

    collider_bias_type <- "slopehunter"
    pruned_harmonised_effects_df <- data.frame(pruned_harmonised_effects)
    slopehunter_result <- SlopeHunter::hunt(
      dat = pruned_harmonised_effects_df,
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

    collider_bias_results <- collider_bias_results %>% add_row(
      METHOD = collider_bias_type,
      P_VAL_THRESHOLD = p_value_threshold,
      BETA = slopehunter_estimated_slope,
      SE = slopehunter_estimated_slope_standard_error,
      HUNTED = hunted,
      ENTROPY = entropy,
      PLEIOTROPIC = pleiotropic
    )

    harmonised_effects <- adjust_gwas_data_from_weights(harmonised_effects, collider_bias_type, slopehunter_estimated_slope, slopehunter_estimated_slope_standard_error)
  }

  if (include_dudbridge) {
    suppressPackageStartupMessages(library(MendelianRandomization))
    print("Starting Dudbridge")

    collider_bias_type <- "dudbridge"

    weighted_dudbridge <- MendelianRandomization::mr_ivw(mr_input(
      bx = pruned_harmonised_effects$BETA.incidence,
      bxse = pruned_harmonised_effects$SE.incidence,
      by = pruned_harmonised_effects$BETA.prognosis,
      byse = pruned_harmonised_effects$SE.prognosis
    ))

    pruned_harmonised_effects$weights <- 1 / pruned_harmonised_effects$SE.prognosis^2

    weighting <- (sum(pruned_harmonised_effects$weights * pruned_harmonised_effects$BETA.incidence^2)) /
      (
        (sum(pruned_harmonised_effects$weights * pruned_harmonised_effects$BETA.incidence^2)) -
          (sum(pruned_harmonised_effects$weights * pruned_harmonised_effects$SE.incidence^2))
      )

    dudbridge_estimated_slope <- weighted_dudbridge$Estimate * weighting
    dudbridge_estimated_slope_standard_error <- weighted_dudbridge$StdError * weighting

    harmonised_effects <- adjust_gwas_data_from_weights(harmonised_effects, collider_bias_type, dudbridge_estimated_slope, dudbridge_estimated_slope_standard_error)

    collider_bias_results <- collider_bias_results %>% add_row(
      METHOD = collider_bias_type,
      P_VAL_THRESHOLD = p_value_threshold,
      BETA = dudbridge_estimated_slope,
      SE = dudbridge_estimated_slope_standard_error,
      HUNTED = NA,
      ENTROPY = NA,
      PLEIOTROPIC = NA
    )
  }

  collider_bias_result_file <- gsub("([\\w/]*)(\\..*)$", "\\1_cb_results\\2", output_file)
  data.table::fwrite(collider_bias_results, collider_bias_result_file)
  # TODO: should we save a file for each collider bias correction?
  data.table::frwite(harmonised_effects, output_file, sep = "\t")
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
adjust_gwas_data_from_weights <- function(gwas, type_of_adjustment, slope, slope_standard_error) {
  type_of_adjustment <- tolower(type_of_adjustment)
  adjusted_beta <- paste0("BETA.", type_of_adjustment)
  adjusted_se <- paste0("SE.", type_of_adjustment)
  adjusted_p <- paste0("P.", type_of_adjustment)

  gwas[[adjusted_beta]] <- (gwas$BETA.prognosis - (slope * gwas$BETA.incidence))

  gwas[[adjusted_se]] <- sqrt(
    (gwas$SE.prognosis * gwas$SE.prognosis) +
      ((slope * slope) * (gwas$SE.incidence * gwas$SE.incidence)) +
      ((gwas$BETA.incidence * gwas$BETA.incidence) * (slope_standard_error * slope_standard_error)) +
      ((gwas$SE.incidence * gwas$SE.incidence) * (slope_standard_error * slope_standard_error))
  )

  gwas[[adjusted_p]] <- pchisq(
    (gwas[[adjusted_beta]] / gwas[[adjusted_se]])^2, 1,
    lower.tail = FALSE
  )

  return(gwas)
}

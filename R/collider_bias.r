collider_bias_type <- list(
  slopehunter = "slopehunter",
  cwls = "cwls",
  mr_ivw = "mr_ivw"
)

collider_bias_results <- data.frame(
  METHOD = character(),
  P_VALUE_THRESHOLD = numeric(),
  SNPS_USED = numeric(),
  BETA = numeric(),
  SE = numeric(),
  PLEIOTROPIC = numeric(),
  ENTROPY = numeric()
)

#' correct_for_collider_bias: do a bunch of collider bias corrections.  Including:
#'  * SlopeHunter
#'  * "Corrected Weighted Least Squares" (CWLS, or Dudbridge) Correction
#'  * MR IVW
#'
#' @return 2 plots: one manhattan plot and one QQ plot (with lambda included)
#' @examples
conduct_collider_bias_analysis <- function(incidence_gwas,
                                           subsequent_gwas,
                                           clumped_snps_file,
                                           collider_bias_results_file,
                                           harmonised_effects_result_file,
                                           slopehunter_adjusted_file,
                                           p_value_thresholds = c(0.1, 0.01, 0.001, 1e-05)) {
  library(dplyr, quietly = TRUE)
  library(SlopeHunter, quietly = TRUE)

  clumped_snps <- data.table::fread(clumped_snps_file)
  incidence <- data.table::fread(incidence_gwas)
  subsequent <- data.table::fread(subsequent_gwas)
  clumped_snps <- map_rsid_list_to_snps(incidence, clumped_snps$SNP)

  print(paste("Found", length(clumped_snps), "in incidence GWAS for clumping"))

  slopehunter_incidence <- SlopeHunter::read_incidence(incidence_gwas,
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

  slopehunter_incidence <- slopehunter_incidence[
    !(is.na(slopehunter_incidence$EA.incidence) | is.na(slopehunter_incidence$OA.incidence)),
  ]
  slopehunter_progression <- SlopeHunter::read_prognosis(subsequent_gwas,
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

  slopehunter_progression <- slopehunter_progression[
    !(is.na(slopehunter_progression$EA.prognosis) | is.na(slopehunter_progression$OA.prognosis)),
  ]

  harmonised_effects <- SlopeHunter::harmonise_effects(
    incidence_dat = slopehunter_incidence,
    prognosis_dat = slopehunter_progression,
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
  vroom::vroom_write(harmonised_effects, harmonised_effects_result_file)

  clumped_snps <- tolower(clumped_snps)

  pruned_harmonised_effects <- subset(harmonised_effects, (SNP %in% clumped_snps) == TRUE)
  pruned_harmonised_effects <- subset(pruned_harmonised_effects, is.na(BETA.prognosis) == FALSE)
  pruned_harmonised_effects <- subset(pruned_harmonised_effects, BETA.prognosis < 10)

  print("Starting SlopeHunter")
  pruned_harmonised_effects_df <- data.frame(pruned_harmonised_effects)

  for (p_value_threshold in p_value_thresholds) {
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

    slopehunter_beta <- slopehunter_result$b
    slopehunter_se <- slopehunter_result$bse
    fit <- slopehunter_result$Fit

    collider_bias_results <- collider_bias_results %>% dplyr::add_row(
      METHOD = collider_bias_type$slopehunter,
      P_VALUE_THRESHOLD = p_value_threshold,
      SNPS_USED = table(fit$clusters)[1],
      BETA = slopehunter_beta,
      SE = slopehunter_se,
      ENTROPY = slopehunter_result$entropy,
      PLEIOTROPIC = table(fit$clusters)[2]
    )
  }

  library(MendelianRandomization, quietly = T)
  print("Starting Dudbridge cwls")

  cwls_coreection <- MendelianRandomization::mr_ivw(MendelianRandomization::mr_input(
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

  cwls_estimated_slope <- cwls_coreection$Estimate * weighting
  cwls_estimated_standard_error <- cwls_coreection$StdError * weighting

  collider_bias_results <- collider_bias_results %>% dplyr::add_row(
    METHOD = collider_bias_type$cwls,
    P_VALUE_THRESHOLD = NA,
    SNPS_USED = length(clumped_snps),
    BETA = cwls_estimated_slope,
    SE = cwls_estimated_standard_error,
    ENTROPY = NA,
    PLEIOTROPIC = NA
  )

  print("Starting TwoSampleMR ivw")
  library(TwoSampleMR, quietly = T)
  for (p_value_threshold in p_value_thresholds) {
    mr_incidence <- TwoSampleMR::read_exposure_data(incidence_gwas,
                                                    snp_col="SNP",
                                                    effect_allele_col = "EA",
                                                    other_allele_col = "OA",
                                                    eaf_col="EAF",
                                                    pval_col="P",
                                                    beta_col="BETA",
                                                    se_col="SE",
                                                    chr_col="CHR",
                                                    pos_col="BP",
                                                    sep="\t"
    )

    mr_incidence <- subset(mr_incidence, pval.exposure <= p_value_threshold) %>%
      subset((SNP %in% clumped_snps) == TRUE) %>%
      subset(duplicated(mr_incidence) == FALSE)
    mr_incidence$exposure <- "incidence"
    mr_incidence$id.exposure <- "incidence"
    mr_incidence$mr_keep.exposure <- TRUE

    mr_progression <- TwoSampleMR::read_outcome_data(file=subsequent_gwas,
                                                     snp_col="SNP",
                                                     effect_allele_col = "EA",
                                                     other_allele_col = "OA",
                                                     eaf_col="EAF",
                                                     pval_col="P",
                                                     beta_col="BETA",
                                                     se_col="SE",
                                                     chr_col="CHR",
                                                     pos_col="BP",
                                                     snps = mr_incidence$SNP,
                                                     sep="\t"
    )

    mr_harmonised <- TwoSampleMR::harmonise_data(mr_incidence, mr_progression)

    mr_results <- TwoSampleMR::mr(mr_harmonised, method_list = c("mr_ivw"))
    collider_bias_results <- collider_bias_results %>% dplyr::add_row(
      METHOD = collider_bias_type$mr_ivw,
      P_VALUE_THRESHOLD = p_value_threshold,
      SNPS_USED = mr_results$nsnp,
      BETA = mr_results$b,
      SE = mr_results$se,
      ENTROPY = NA,
      PLEIOTROPIC = NA
    )
  }

  vroom::vroom_write(collider_bias_results, collider_bias_results_file)

  slopehunter_default_result <- subset(collider_bias_results, METHOD == collider_bias_type$slopehunter & P_VALUE_THRESHOLD == 0.001)
  harmonised_effects <- adjust_gwas_data_from_weights_and_save(subsequent,harmonised_effects,
                                                               slopehunter_default_result$METHOD,
                                                               slopehunter_default_result$BETA,
                                                               slopehunter_default_result$SE,
                                                               slopehunter_adjusted_file
  )
}

#' adjust_gwas_data_from_weights: Apply a slope correction weight (and SE) to an existing GWAS and save result to disk.
#'   This can be used in conjuction with weights calculated to account for collider bias.
#'   Currently used to work on a GWAS result of Slopehunter.
#'
#' @param gwas: gwas that the new collider bias correction will be saved to.  Swaps out BETA, SE, and P
#' @param harmonised_effects: a dataframe that includes BETA.incidence, BETA.prognosis, SE.incidence, SE.prognosis
#' @param collider_bias_type: string name of adjustment (eg. slopehunter)
#' @param beta: number, slope of the correction
#' @param se: number, SE of the corrected slope
#' @param output_file: name of file that the corrected GWAS will be saved to
#'
adjust_gwas_data_from_weights_and_save <- function(gwas,
                                                   harmonised_effects,
                                                   collider_bias_type,
                                                   beta,
                                                   se,
                                                   output_file) {
  harmonised_effects$SNP <- toupper(harmonised_effects$SNP)
  adjusted_beta <- paste0("BETA.", collider_bias_type)
  adjusted_se <- paste0("SE.", collider_bias_type)
  adjusted_p <- paste0("P.", collider_bias_type)

  harmonised_effects[[adjusted_beta]] <- (harmonised_effects$BETA.prognosis - (beta * harmonised_effects$BETA.incidence))

  harmonised_effects[[adjusted_se]] <- sqrt(
    (harmonised_effects$SE.prognosis * harmonised_effects$SE.prognosis) +
      ((beta * beta) * (harmonised_effects$SE.incidence * harmonised_effects$SE.incidence)) +
      ((harmonised_effects$BETA.incidence * harmonised_effects$BETA.incidence) * (se * se)) +
      ((harmonised_effects$SE.incidence * harmonised_effects$SE.incidence) * (se * se))
  )

  harmonised_effects[[adjusted_p]] <- pchisq(
    (harmonised_effects[[adjusted_beta]] / harmonised_effects[[adjusted_se]])^2, 1,
    lower.tail = FALSE
  )

  gwas <- merge(gwas, subset(harmonised_effects, select = c("SNP", adjusted_beta, adjusted_se, adjusted_p)), by="SNP")
  gwas <- subset(gwas, select = -c(BETA, SE, P)) %>%
    rename(BETA = adjusted_beta, SE = adjusted_se, P = adjusted_p)

  vroom::vroom_write(gwas, output_file)
}
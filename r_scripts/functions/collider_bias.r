collider_bias_type <- list(
  slopehunter = "slopehunter",
  dudbridge = "dudbridge",
  ivw = "ivw"
)

collider_bias_results <- data.frame(
  METHOD = character(),
  P_VAL_THRESHOLD = numeric(),
  BETA = numeric(),
  SE = numeric(),
  HUNTED = numeric(),
  PLEIOTROPIC = numeric(),
  ENTROPY = numeric()
)

#' correct_for_collider_bias: do a bunch of collider bias corrections.  Including:
#'  * SlopeHunter
#'  * Dudbridge Correction
#'  * MR IVW
#'
#' @return 2 plots: one manhattan plot and one QQ plot (with lambda included)
#' @examples
correct_for_collider_bias <- function(incidence_gwas,
                                      subsequent_gwas,
                                      clumped_snps_file,
                                      collider_bias_results_file,
                                      slopehunter_adjusted_file,
                                      dudbridge_adjusted_file,
                                      #TODO: maybe change p_value_theshold to list?
                                      p_value_threshold = 0.001,
                                      include_slopehunter = T,
                                      include_dudbridge = T) {
  library(dplyr, quietly = TRUE)
  library(SlopeHunter, quietly = TRUE)

  clumped_snps <- data.table::fread(clumped_snps_file)
  incidence <- data.table::fread(incidence_gwas)
  subsequent <- data.table::fread(subsequent_gwas)
  clumped_snps <- map_rsid_list_to_snps(incidence, clumped_snps$SNP)
  print(paste("Found", length(clumped_snps), "in incidence GWAS for clumping"))

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
    other_allele_col = "OA",
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
  clumped_snps <- tolower(clumped_snps)

  pruned_harmonised_effects <- subset(harmonised_effects, (SNP %in% clumped_snps) == TRUE)
  pruned_harmonised_effects <- subset(pruned_harmonised_effects, is.na(BETA.prognosis) == FALSE)
  pruned_harmonised_effects <- subset(pruned_harmonised_effects, BETA.prognosis < 10)

  if (include_slopehunter) {
    print("Starting SlopeHunter")

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

    slopehunter_beta <- slopehunter_result$b
    slopehunter_se <- slopehunter_result$bse
    fit <- slopehunter_result$Fit

    collider_bias_results <- collider_bias_results %>% dplyr::add_row(
      METHOD = collider_bias_type$slopehunter,
      P_VAL_THRESHOLD = p_value_threshold,
      BETA = slopehunter_beta,
      SE = slopehunter_se,
      HUNTED = table(fit$clusters)[1],
      ENTROPY = slopehunter_result$entropy,
      PLEIOTROPIC = table(fit$clusters)[2]
    )

    harmonised_effects <- adjust_gwas_data_from_weights(harmonised_effects, collider_bias_type$slopehunter, slopehunter_beta, slopehunter_se)
    save_subsequent_adjusted(collider_bias_type$slopehunter, subsequent, harmonised_effects, slopehunter_adjusted_file)
  }

  if (include_dudbridge) {
    suppressPackageStartupMessages(library(MendelianRandomization))
    print("Starting Dudbridge")

    weighted_dudbridge <- MendelianRandomization::mr_ivw(MendelianRandomization::mr_input(
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

    harmonised_effects <- adjust_gwas_data_from_weights(harmonised_effects, collider_bias_type$dudbridge, dudbridge_estimated_slope, dudbridge_estimated_slope_standard_error)
    save_subsequent_adjusted(collider_bias_type$dudbridge, subsequent, harmonised_effects, dudbridge_adjusted_file)

    collider_bias_results <- collider_bias_results %>% dplyr::add_row(
      METHOD = collider_bias_type$dudbridge,
      P_VAL_THRESHOLD = p_value_threshold,
      BETA = dudbridge_estimated_slope,
      SE = dudbridge_estimated_slope_standard_error,
      HUNTED = NA,
      ENTROPY = NA,
      PLEIOTROPIC = NA
    )
  }
  # if (include_ivw) {
  #   library(TwoSampleMR, quietly = T)
  #   incidence_mr <- TwoSampleMR::read_exposure_data(incidence_gwas,
  #                                    sep="\t",
  #                                    snp_col="SNP",
  #                                    effect_allele_col = "ALLELE1",
  #                                    other_allele_col = "ALLELE0",
  #                                    eaf_col="A1FREQ",
  #                                    pval_col="P_BOLT_LMM_INF",
  #                                    beta_col="BETA",
  #                                    se_col="SE",
  #                                    chr_col="CHR",
  #                                    pos_col="BP"
  #   )
  #   incidence_mr<-subset(incidence_mr, pval.exposure<=5e-8)
  #   incidence_mr<-subset(incidence_mr, (SNP %in% clumped_snps) == TRUE)
  #   incidence_mr$exposure <- "incidence"
  #   incidence_mr$id.exposure <- "incidence"
  #   incidence_mr <- subset(incidence_mr, duplicated(incidence_mr) == FALSE)
  #   incidence_mr$mr_keep.exposure <- TRUE
  #
  #   progression<-read_outcome_data(file=subsequent_gwas,
  #                                  snp_col="SNP",
  #                                  snps = incidence_mr$SNP,
  #                                  effect_allele_col = "EA",
  #                                  other_allele_col = "OA",
  #                                  pval_col="P",
  #                                  beta_col="BETA",
  #                                  se_col="SE",
  #                                  chr_col="CHR",
  #                                  pos_col="BP",
  #                                  sep="\t"
  #   )
  #
  #   progression<-subset(progression, select = -c(eaf.outcome))
  #   progression<-merge(progression, af_plink, by="SNP")
  #   dat<-harmonise_data(incidence_mr, progression)
  #   dat$beta.outcome<-log(dat$beta.outcome)
  #
  #   mr_res<-TwoSampleMR::mr(dat, method_list = c("mr_ivw","mr_egger_regression", "mr_weighted_median","mr_weighted_mode"))
  #   write.table(mr_res, file=paste(filenames[k],"_MR_results.txt",sep=""), row.names=F, quote=F)
  #
  #   cf.mr<-mr_res[1,7]
  #   cf.se.mr<-mr_res[1,8]
  #
  # }
  data.table::fwrite(collider_bias_results, collider_bias_results_file, sep="\t")
}

#' adjust_gwas_data_from_weights: Apply a slope correction weight (and SE) to an existing GWAS
#'   This can be used in conjuction with weights calculated to account for collider bias.
#'   Currently used to work on a GWAS result of Slopehunter.
#'
#' @param gwas: a dataframe that includes BETA.incidence, BETA.prognosis, SE.incidence, SE.prognosis
#' @param collider_bias_type: string name of adjustment (eg. SLOPEHUNTER)
#' @param slope: number, slope of the correction
#' @param slope_standard_error: number, SE of the corrected slope
#' @return gwas with 3 additional columns BETA_{name}, SE_{name}, and P_{name}
adjust_gwas_data_from_weights <- function(gwas, collider_bias_type, slope, slope_standard_error) {
  adjusted_beta <- paste0("BETA.", collider_bias_type)
  adjusted_se <- paste0("SE.", collider_bias_type)
  adjusted_p <- paste0("P.", collider_bias_type)

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

save_subsequent_adjusted <- function(collider_bias_type, gwas, harmonised_effects, output_file) {
  gwas$P <- harmonised_effects[[paste0("P.", collider_bias_type)]]
  gwas$SE <- harmonised_effects[[paste0("SE.", collider_bias_type)]]
  gwas$BETA <- harmonised_effects[[paste0("BETA.", collider_bias_type)]]

  data.table::fwrite(gwas, output_file, sep = "\t")
}
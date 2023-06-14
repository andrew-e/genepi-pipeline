
read_and_format <-function(file_gwas, data_version, custom_file_columns=NA){
  
  # different versions of data formats to read in
  if (data_version == "ieu_ukb_gwas_pipeline"){
    # data produced by IEU GWAS pipeline
    out <-vroom(file_gwas,
                col_select = c("SNP","BETA","SE","ALLELE1","ALLELE0","A1FREQ","P_BOLT_LMM_INF")) %>% 
      format_data(., type="outcome",
                  snp_col = "SNP",
                  beta_col = "BETA",
                  se_col = "SE",
                  effect_allele_col = "ALLELE1",
                  other_allele_col = "ALLELE0",
                  eaf_col = "A1FREQ",
                  pval_col = "P_BOLT_LMM_INF")
  
    

  }else if (data_version == "ukb_neale"){
    # UKB data produced by Neale lab
    print("reading variants")
    variants <- vroom(paste0("meta/neale_lab_variants_rsid_only.tsv")) # this is a col subset of variants.tsv.bgz from Neale lab
    print("reading gwas")
    gwas <-vroom(file_gwas, col_select = c("variant","minor_allele","minor_AF","beta","se","pval"))
    
    print("joining")
    if (nrow(gwas) == nrow(variants)){
      merged <- left_join(gwas, variants, by = c("variant"="variant", "minor_allele"="alt")) %>% 
        select(-variant) %>% select(SNP=rsid, everything())
    }
    
    print("formatting")
    out <-format_data(merged, type="outcome",
                      snp_col = "SNP",
                      beta_col = "beta",
                      se_col = "se",
                      effect_allele_col = "minor_allele",
                      other_allele_col = "ref",
                      eaf_col = "minor_AF",
                      pval_col = "pval")
    
    
  } else if (data_version == 'gwas_cat'){
    # harmonised data fromat from GWAS catalog
    out <-vroom(file_gwas,
                col_select = c("variant_id",'beta', "standard_error","effect_allele","other_allele","effect_allele_frequency","p_value")) %>% 
      format_data(., type="outcome",
                  snp_col = "variant_id",
                  beta_col = "beta",
                  se_col = "standard_error",
                  effect_allele_col = "effect_allele",
                  other_allele_col = "other_allele",
                  eaf_col = "effect_allele_frequency",
                  pval_col = "p_value")
    
    
  } else if (data_version == "custom") {
    
    
    custom_file_columns <- custom_file_columns %>% select_if(~all(!is.na(.)))
    stopifnot(dim(custom_file_columns)[2]>0)
    
    
    # data produced by IEU GWAS pipeline: pval col P_BOLT_LMM
    out <-vroom(file_gwas,
                col_select = c(custom_file_columns[1,])) %>% 
      format_data(., 
                  phenotype_col= custom_file_columns$trait,
                  snp_col = custom_file_columns$custom_SNP,
                  beta_col = custom_file_columns$custom_beta,
                  se_col = custom_file_columns$custom_se,
                  effect_allele_col = custom_file_columns$custom_effect_allele,
                  other_allele_col = custom_file_columns$custom_other_allele,
                  eaf_col = custom_file_columns$custom_eaf,
                  pval_col = custom_file_columns$custom_pval,
                  samplesize_col = custom_file_columns$custom_samplesize,
                  ncase_col = custom_file_columns$custom_ncase ,
                  ncontrol_col = custom_file_columns$custom_ncontrol,
                  chr_col = custom_file_columns$custom_chr,
                  pos_col = custom_file_columns$custom_pos,
                  z_col = custom_file_columns$custom_z,
                  info_col = custom_file_columns$custom_info, 
                  units_col = custom_file_columns$custom_units,
                  gene_col = custom_file_columns$custom_gene
                  ) %>% 
      #TODO drop .outcome
      
    
    
    
    
    
    
  }
  return(out)
  
}

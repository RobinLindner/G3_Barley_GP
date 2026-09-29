source("0_utils.R")
## ---- Computing BLUEs - Regular traits ----
file.remove(BLUE_variance_components_path)
file.remove(BLUE_path)
if(!(file.exists(BLUE_variance_components_path) & file.exists(BLUE_path))) {
  
  pheno_long = read.csv(phenotype_nonHSR_long_file)
  t_d_comb = pheno_long %>%
    dplyr::select(Trait,DAT) %>%
    distinct(.keep_all = T)
  
  Var_comp_frame = data.frame(DAT=NA,Trait=NA,Expvar=NA,GxEvar=NA,Repvar=NA,evar=NA)
  BLUE_frame = data.frame(Genotype=NA,Trait=NA,DAT=NA,Value=NA)
  for(i in 1:nrow(t_d_comb)){
    print(paste0("Computing BLUEs for trait: ",t_d_comb$Trait[i]," at ",t_d_comb$DAT[i]," DAT."))
    unf_selection <- pheno_long %>%
      filter(DAT==t_d_comb$DAT[i] & Trait == t_d_comb$Trait[i]) 
    
    selection <- unf_selection %>%
      filter(!findExtreme(Value))
    
    diff = nrow(selection) - nrow(unf_selection)
    if(diff!=0){
      print(paste0("Filtered out: ",diff," extreme observations."))
    }
    
    # remove measurement instances with less than 100 measured plants
    if(nrow(selection)<100){next}
    
    # For these groups the pots were measured, removing the estimated replicate plant effect (var(rep)=0).
    if(unique(selection$Group)=="SW"| unique(selection$Group)=="RGB"){
      model_formula = as.formula(Value ~ 0 + Genotype + (1|ExpID) + (1|Genotype:ExpID))
      
    }
    # For these cases, only a single time point was measured, effectively removing estimable environmental variance (var(E)=0)
    else if(var(selection$ExpID)==0){
      model_formula = as.formula(Value ~ 0 + Genotype + (1|Genotype:Rep_ID))
      
    }
    # in other cases, the full model as described in the main text was fitted.
    else{
      model_formula = as.formula(Value ~ 0 + Genotype + (1|Genotype:Rep_ID) + (1|ExpID) + (1|Genotype:ExpID))
    }
    BLUE_model = lmer(model_formula,
                      data = selection,
                      control=lmerControl(optimizer ="bobyqa",
                                          check.conv.singular = .makeCC(action = "message",  tol = 1e-4))) # define model
    
    random_effects <- ranef(BLUE_model)
    fixed_effects <- fixef(BLUE_model)
    BLUEs = fixed_effects
    Genotypes = gsub("Genotype","",names(BLUEs))
    
    var.comp=VarCorr(BLUE_model)
    if(is.null(var.comp$ExpID)){
      Expvar = 0
    }else{
      Expvar= var.comp$ExpID 
    }
    if(is.null(var.comp$`Genotype:ExpID`)){
      GxEvar = 0
    }else{
      GxEvar= var.comp$`Genotype:ExpID`  
    }
    if(is.null(var.comp$`Genotype:Rep_ID`)){
      Repvar = 0
    }else{
      Repvar= var.comp$`Genotype:Rep_ID` 
    }
    evar= attr(var.comp,'sc')^2
    
    
    cur_frame = data.frame(Genotype=Genotypes,
                           Trait=t_d_comb$Trait[i],
                           DAT=t_d_comb$DAT[i],
                           Value=BLUEs)
    
    rownames(cur_frame) = NULL
    
    Var_comp_frame[i,] = c(
      t_d_comb$DAT[i],
      t_d_comb$Trait[i],
      Expvar,
      GxEvar,
      Repvar,
      evar)
    if(i==1){
      BLUE_frame = cur_frame
    }else{
      BLUE_frame = rbind(BLUE_frame,cur_frame)
    }
    
  }
  
  Var_comp_frame <- Var_comp_frame%>%
    filter(!is.na(DAT))
  
  # Write out the variance components and BLUPs for each measurement instance (trait-time point pair)
  write.csv(Var_comp_frame, BLUE_variance_components_path,row.names = F)
  write.csv(BLUE_frame,BLUE_path,row.names=F)
}else{
  Var_comp_frame = read.csv(BLUE_variance_components_path) 
  BLUP_frame = read.csv(BLUE_path) 
}
source("0_utils.R")
file.remove(HSR_BLUE_variance_components_path)
file.remove(HSR_BLUE_path)
if(!(file.exists(HSR_BLUE_variance_components_path) & 
     file.exists(HSR_BLUE_path))) {
  full_spectrum = read.csv(phenotype_HSR_file)
  
  long_HS <- full_spectrum %>%
    mutate(ExpID = substring(ID,1,1),Rep_ID = substring(ID,12,12)) %>%
    pivot_longer(cols=colnames(full_spectrum)[5:ncol(full_spectrum)],names_to = "Trait",values_to="Value") %>%
    filter(!is.na(Value)) 
  
  t_d_comb_HS = long_HS %>%
    dplyr::select(Trait,DAT) %>%
    distinct(.keep_all = T)
  
  Var_comp_frame_HS = data.frame(DAT=NA,Trait=NA,Expvar=NA,GxEvar=NA,Repvar=NA,evar=NA)
  BLUE_frame_HS = data.frame(Genotype=NA,Trait=NA,DAT=NA,Value=NA)
  
  for(i in 1:nrow(t_d_comb_HS)){
    print(paste0("Computing BLUEs for trait: ",t_d_comb_HS$Trait[i],
                 " at ",t_d_comb_HS$DAT[i]," DAT."))
    
    unf_selection <- long_HS %>%
      filter(DAT==t_d_comb_HS$DAT[i] & Trait == t_d_comb_HS$Trait[i])
    
    selection <- unf_selection %>%
      filter(!findExtreme(Value))
    
    if(var(selection$ExpID)==0){
      model_formula = as.formula(Value ~ 0 + Genotype + (1|Genotype:Rep_ID))
      
    }else{
      model_formula = as.formula(Value ~ 0 + Genotype + (1|Genotype:Rep_ID) + (1|ExpID) + (1|Genotype:ExpID))
    }
    
    BLUE_model=lmer(model_formula,
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
                           Trait=t_d_comb_HS$Trait[i],
                           DAT=t_d_comb_HS$DAT[i],
                           Value=BLUEs)
    
    rownames(cur_frame) = NULL
    
    Var_comp_frame_HS[i,] = c(
      t_d_comb_HS$DAT[i],
      t_d_comb_HS$Trait[i],
      Expvar,
      GxEvar,
      Repvar,
      evar)
    if(i==1){
      BLUE_frame_HS = cur_frame
    }else{
      BLUE_frame_HS = rbind(BLUE_frame_HS,cur_frame)
    }
    
  }
  
  Var_comp_frame_HS <- Var_comp_frame_HS %>%
    mutate(WL=as.numeric(sub("VNIR_Plant_Spectrum_X","",Trait)))
  
  Var_comp_frame_HS$Trait = t_d_comb_HS$Trait
  Var_comp_frame_HS$DAT = t_d_comb_HS$DAT
  Var_comp_frame_HS$DAT = as.numeric(Var_comp_frame_HS$DAT)
  Var_comp_frame_HS$Expvar = as.numeric(Var_comp_frame_HS$Expvar)
  Var_comp_frame_HS$GxEvar = as.numeric(Var_comp_frame_HS$GxEvar)
  Var_comp_frame_HS$Repvar = as.numeric(Var_comp_frame_HS$Repvar)
  Var_comp_frame_HS$evar = as.numeric(Var_comp_frame_HS$evar)
  
  
  
  write.csv(Var_comp_frame_HS,HSR_BLUE_variance_components_path,row.names = F)
  write.csv(BLUE_frame_HS,HSR_BLUE_path,row.names = F)
  
}else{
  Var_comp_frame_HS = read.csv(HSR_BLUE_variance_components_path) 
  BLUE_frame_HS = read.csv(HSR_BLUE_path) 
}
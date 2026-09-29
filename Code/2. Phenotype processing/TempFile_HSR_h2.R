source("0_utils.R")
K = read.csv(GRM_path,row.names = 1)

full_spectrum = read.csv(phenotype_HSR_file)
long_HS <- full_spectrum %>%
  mutate(ExpID = substring(ID,1,1),Rep_ID = substring(ID,12,12)) %>%
  pivot_longer(cols=colnames(full_spectrum)[5:ncol(full_spectrum)],names_to = "Trait",values_to="Value") %>%
  filter(!is.na(Value)) 
t_d_comb = long_HS %>%
  dplyr::select(Trait,DAT) %>%
  distinct(.keep_all = T)

#BLUE_frame_HS = read.csv(HSR_BLUE_path)

for(i in 1:nrow(t_d_comb)){
  print(paste0("Computing h^2 for trait: ",t_d_comb$Trait[i]," at ",t_d_comb$DAT[i]," DAT."))
  unf_selection <- long_HS %>%
    filter(DAT==t_d_comb$DAT[i] & Trait == t_d_comb$Trait[i] & Genotype %in% rownames(K)) 
  
  # We filter for measurement errors by a threshold of  Q1(x)-10*sd(x)< x <Q3(x) + 10*sd(x)
  selection <- unf_selection %>%
    filter(!findExtreme(Value))
  
  diff = nrow(selection) - nrow(unf_selection)
  if(diff!=0){
    print(paste0("Filtered out: ",diff," extreme observations."))
  }
  
  # remove measurement instances with less than 100 measured plants
  if(nrow(selection)<100){next}
  
  # For these groups the pots were measured, removing the estimated replicate plant effect (var(rep)=0).
  if(var(selection$Rep_ID)==0){
    model_formula = as.formula(Value ~ (1|Genotype) + (1|ExpID) + (1|Genotype:ExpID))
  }
  # For these cases, only a single time point was measured, effectively removing estimable environmental variance (var(E)=0)
  else if(var(selection$ExpID)==0){
    model_formula = as.formula(Value ~ (1|Genotype)  + (1|Genotype:Rep_ID))
  }
  # in other cases, the full model as described in the main text was fitted.
  else{
    model_formula = as.formula(Value ~ (1|Genotype)  + (1|Genotype:Rep_ID) + (1|ExpID) + (1|Genotype:ExpID))
  }
  BLUP_model = relmatLmer(model_formula,
                          data = selection,
                          relmat = list(Genotype=as.matrix(K)),
                          control=lmerControl(optimizer ="bobyqa",
                                              check.conv.singular = .makeCC(action = "message",  tol = 1e-4))) # define model
  
  random_effects <- ranef(BLUP_model)
  
  var.comp=VarCorr(BLUP_model)
  
  Gvar = var.comp$Genotype
  
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
  
  
  
  h2 = Gvar / (Gvar+Expvar+Repvar+GxEvar+evar)
  
  cur_df = data.frame(Trait = t_d_comb$Trait[i],
                      DAT = t_d_comb$DAT[i],
                      Gvar = Gvar,
                      Expvar = Expvar,
                      Repvar = Repvar,
                      GxEvar = GxEvar,
                      evar = evar,
                      h2 = h2)
  
  colnames(cur_df) = c("Trait","DAT","Gvar","Expvar","Repvar","GxEvar","evar","h2")
  rownames(cur_df) = NULL
  
  if(i==1){
    h2_df_HS = cur_df
  }else{
    h2_df_HS = rbind(h2_df_HS,cur_df)
  }
}

write.csv(h2_df_HS,HSR_h2_path)

source("0_utils.R")
## ---- Estimating narrow-sense heritability ---- 

pheno_long = read.csv(phenotype_nonHSR_long_file)

t_d_comb = pheno_long %>%
  dplyr::select(Trait,DAT) %>%
  distinct(.keep_all = T)

K = read.csv(GRM_path,row.names = 1)
missingTnD = c()
i=30
for(i in 1:nrow(t_d_comb)){
  print(paste0("Computing h^2 for trait: ",t_d_comb$Trait[i]," at ",t_d_comb$DAT[i]," DAT."))
  unf_selection <- pheno_long %>%
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
  if(unique(selection$Group)=="SW"| unique(selection$Group)=="RGB"){
    model_formula = as.formula(Value ~ (1|Genotype) + (1|ExpID) + (1|Genotype:ExpID))
    
  }else if(var(selection$ExpID)==0){ # For these cases, only a single time point was measured, effectively removing estimable environmental variance (var(E)=0)
    model_formula = as.formula(Value ~ (1|Genotype)  + (1|Genotype:Rep_ID))
  }else{# in other cases, the full model as described in the main text was fitted.
    model_formula = as.formula(Value ~ (1|Genotype)  + (1|Genotype:Rep_ID) + (1|ExpID) + (1|Genotype:ExpID))
  }
  BLUP_model = relmatLmer(model_formula,
                          data = selection,
                          relmat = list(Genotype=as.matrix(K)),
                          control=lmerControl(optimizer ="bobyqa",
                                              check.conv.singular = .makeCC(action = "message",  tol = 1e-4))) # define model
  
  
  harmonic_mean <- function(x, na.rm = FALSE) {
    if (na.rm) x <- x[!is.na(x)]
    if (any(x == 0)) stop("Harmonic mean is undefined when any value is zero")
    n <- length(x)
    n / sum(1 / x)
  }
  n_E = BLUP_model@frame %>% 
    group_by(Genotype) %>%
    summarize(count = length(unique(ExpID)))
  harmonic_mean(n_E$count,na.rm=T)

  n_Rep = BLUP_model@frame %>% 
  group_by(Genotype) %>%
  summarize(count = length(unique(Rep_ID)))  
  
  
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
  
  p_mean = mean(selection$Value,na.rm=T)
  
  h2 = Gvar / (Gvar+Expvar+Repvar+GxEvar+evar)
  
  cur_df = data.frame(Trait = t_d_comb$Trait[i],
                      DAT = t_d_comb$DAT[i],
                      Gvar = Gvar,
                      Expvar = Expvar,
                      Repvar = Repvar,
                      GxEvar = GxEvar,
                      evar = evar,
                      h2 = h2,
                      pheno_mean = p_mean)
  colnames(cur_df) = c("Trait","DAT","Gvar","Expvar","Repvar","GxEvar","evar","h2","pheno_mean")
  rownames(cur_df) = NULL
  
  if(i==1){
    h2_df = cur_df
  }else{
    h2_df = rbind(h2_df,cur_df)
  }
}

colnames(h2_df) = c("Trait","DAT","Gvar","Expvar","Repvar","GxEvar","evar","h2","pheno_mean")
rownames(h2_df) = NULL

write.csv(h2_df,nonHSR_h2_path)

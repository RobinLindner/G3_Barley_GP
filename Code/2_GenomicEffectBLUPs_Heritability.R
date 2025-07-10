## ---------------------------
##
## Script name: 2_GenomicEffectBLUPs_Heritability
##
## Purpose of script: Contains analyses conducted for named aspects of the G3
##                    Barley time-series GWAS and GP publication.
##
## Author: M.Sc. Robin Lindner
##
## Date Created: 2025-07-03
##
## Copyright (c) Robin Lindner, 2025
## Email: robin.lindner@uni-potsdam.de
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

source("0_utils.R")

## ---------------------------

## ---- Measurement instances ---- 
# Mapping of traits to groups
trait_groups = read.csv(trait_groups_file)

phenotype_nonHSR = read.csv(phenotype_nonHSR_file)

# Remove environmental parameters (last 16 columns)
phenotype_nonHSR_noEvn = phenotype_nonHSR[,1:125]

# Create long data frame, converting Traits to key-value variables 
pheno_long <- phenotype_nonHSR_noEvn %>%
  mutate(ExpID = substring(ID,1,1),Rep_ID = substring(ID,12,12)) %>%
  pivot_longer(cols=colnames(phenotype_nonHSR_noEvn)[5:125],names_to = "Trait",values_to="Value") %>%
  filter(!is.na(Value)) %>%
  mutate(Group=trait_groups$Group[match(Trait,trait_groups$Trait)])

write.csv(pheno_long,file = phenotype_nonHSR_long_file,row.names = F)

# Find trait and time point combinations with at least one observation.
t_d_comb = pheno_long %>%
  dplyr::select(Trait,DAT) %>%
  distinct(.keep_all = T)

# count the number of observations/measurements for each trait-time point pair
measurement_summary = data.frame(Trait=NA,DAT=0,nobs=0,ngeno=0)
for(i in 1:nrow(t_d_comb)){
  selection <- pheno_long %>%
    filter(DAT==t_d_comb$DAT[i] & Trait == t_d_comb$Trait[i])
  measurement_summary[i,1] = t_d_comb$Trait[i]
  measurement_summary[i,2] = as.numeric(t_d_comb$DAT[i])
  measurement_summary[i,3] = nrow(selection)
  measurement_summary[i,4] = length(unique(selection$Genotype))
}


# Find measurements that have a sample size that is too small
sum_geno <-measurement_summary %>%
  group_by(ngeno)%>%
  summarise(count=n()) %>%
  mutate(cs = cumsum(count))

sum_obs <-measurement_summary %>%
  group_by(nobs)%>%
  summarise(count=n()) %>%
  mutate(cs = cumsum(count))

summary(measurement_summary)

# 2 measurements with only 22 measured plants and 18 genotypes are subsequently removed

## ---- Computing BLUPs - Regular traits ----
if(!(file.exists(BLUP_variance_components_path) & file.exists(BLUP_path))) {
Var_comp_frame = data.frame(DAT=NA,Trait=NA,Intercept=NA,Gvar=NA,Expvar=NA,GxEvar=NA,Repvar=NA,evar=NA,H2=NA)
BLUP_frame = data.frame(Genotype=NA,Trait=NA,DAT=NA,Value=NA)
for(i in 1:nrow(t_d_comb)){
  print(paste0("Computing BLUPs for trait: ",t_d_comb$Trait[i]," at ",t_d_comb$DAT[i]," DAT."))
  selection <- pheno_long %>%
    filter(DAT==t_d_comb$DAT[i] & Trait == t_d_comb$Trait[i])
  
  # remove measurement instances with less than 100 measured plants
  if(nrow(selection)<100){next}
  
  # For these groups the pots were measured, removing the estimated replicate plant effect (var(rep)=0).
  if(unique(selection$Group)=="SW"| unique(selection$Group)=="RGB"){
    model=lmer(Value ~ 1 + (1|Genotype) + (1|ExpID) + (1|Genotype:ExpID),
               data = selection,
               control=lmerControl(optimizer ="bobyqa",
                                   check.conv.singular = .makeCC(action = "message",  tol = 1e-4))) # define model
    random_effects <- ranef(model)
    BLUPs = random_effects$Genotype
    Genotypes = rownames(BLUPs)
    
    var.comp=VarCorr(model)
    Inter = fixef(model)
    Gvar =var.comp$Genotype 
    Expvar= var.comp$ExpID 
    GxEvar= var.comp$`Genotype:ExpID` 
    Repvar= 0
    evar= attr(var.comp,'sc')^2
    H2= Gvar/(Gvar+Expvar+GxEvar+Repvar+evar)
  }
  # For these cases, only a single time point was measured, effectively removing estimable environmental variance (var(E)=0)
  else if(var(selection$ExpID)==0){
    model=lmer(Value ~ 1 + (1|Genotype/Rep_ID),
               data = selection,
               control=lmerControl(optimizer ="bobyqa",
                                   check.conv.singular = .makeCC(action = "message",  tol = 1e-4))) # define model
    random_effects <- ranef(model)
    BLUPs = random_effects$Genotype
    Genotypes = rownames(BLUPs)
    
    var.comp=VarCorr(model)
    Inter = fixef(model)
    Gvar =var.comp$Genotype 
    Expvar= 0
    GxEvar= 0
    Repvar= var.comp$`Rep_ID:Genotype` 
    evar= attr(var.comp,'sc')^2
    H2= Gvar/(Gvar+Expvar+GxEvar+Repvar+evar)
  }
  # in other cases, the full model as described in the main text was fitted.
  else{
    model=lmer(Value ~ 1 + (1|Genotype/Rep_ID) + (1|ExpID) + (1|Genotype:ExpID),
               data = selection,
               control=lmerControl(optimizer ="bobyqa",
                                   check.conv.singular = .makeCC(action = "message",  tol = 1e-4))) # define model
    random_effects <- ranef(model)
    BLUPs = random_effects$Genotype
    Genotypes = rownames(BLUPs)
    
    var.comp=VarCorr(model)
    Inter = fixef(model)
    Gvar =var.comp$Genotype 
    Expvar= var.comp$ExpID 
    GxEvar= var.comp$`Genotype:ExpID` 
    Repvar= var.comp$`Rep_ID:Genotype` 
    evar= attr(var.comp,'sc')^2
    H2= Gvar/(Gvar+Expvar+GxEvar+Repvar+evar)
  }
  
  cur_frame = data.frame(Genotype=Genotypes,
                         Trait=t_d_comb$Trait[i],
                         DAT=t_d_comb$DAT[i],
                         Value=BLUPs$`(Intercept)`)
  
  Var_comp_frame[i,] = c(
    t_d_comb$DAT[i],
    t_d_comb$Trait[i],
    Inter,
    Gvar,
    Expvar,
    GxEvar,
    Repvar,
    evar,
    H2)
  if(i==1){
    BLUP_frame = cur_frame
  }else{
    BLUP_frame = rbind(BLUP_frame,cur_frame)
  }
  
}

Var_comp_frame <- Var_comp_frame%>%
  filter(!is.na(DAT))

# Write out the variance components and BLUPs for each measurement instance (trait-time point pair)
write.csv(Var_comp_frame, BLUP_variance_components_path,row.names = F)
write.csv(BLUP_frame,BLUP_path,row.names=F)
}else{
  Var_comp_frame = read.csv(BLUP_variance_components_path) 
  BLUP_frame = read.csv(BLUP_path) 
}



## ---- Computing BLUPs - HSR data ----
if(!(file.exists(HSR_BLUP_variance_components_path) & 
     file.exists(HSR_BLUP_path))) {
  full_spectrum = read.csv(phenotype_HSR_file)
  
  long_HS <- full_spectrum %>%
    mutate(ExpID = substring(ID,1,1),Rep_ID = substring(ID,12,12)) %>%
    pivot_longer(cols=colnames(full_spectrum)[5:ncol(full_spectrum)],names_to = "Trait",values_to="Value") %>%
    filter(!is.na(Value)) 
  
  t_d_comb_HS = long_HS %>%
    dplyr::select(Trait,DAT) %>%
    distinct(.keep_all = T)
  
  Var_comp_frame_HS = data.frame(DAT=NA,Trait=NA,Intercept=NA,Gvar=NA,Expvar=NA,GxEvar=NA,Repvar=NA,evar=NA,H2=NA)
  BLUP_frame_HS = data.frame(Genotype=NA,Trait=NA,DAT=NA,Value=NA)
  for(i in 1:nrow(t_d_comb_HS)){
    print(paste0("Computing BLUPs for trait: ",t_d_comb_HS$Trait[i],
                 " at ",t_d_comb_HS$DAT[i]," DAT."))
    selection <- long_HS %>%
      filter(DAT==t_d_comb_HS$DAT[i] & Trait == t_d_comb_HS$Trait[i])
    
    if(var(selection$ExpID)==0){
      model=lmer(Value ~ 1 + (1|Genotype/Rep_ID),
                 data = selection,
                 control=lmerControl(optimizer ="bobyqa",
                                     check.conv.singular = .makeCC(action = "message",  tol = 1e-4))) # define model
      random_effects <- ranef(model)
      BLUPs = random_effects$Genotype
      Genotypes = rownames(BLUPs)
      
      var.comp=VarCorr(model)
      Inter = fixef(model)
      Gvar =var.comp$Genotype 
      Expvar= 0
      GxEvar= 0
      Repvar= var.comp$`Rep_ID:Genotype` 
      evar= attr(var.comp,'sc')^2
      H2= Gvar/(Gvar+Expvar+GxEvar+Repvar+evar)
    }else{
      model=lmer(Value ~ 1 + (1|Genotype/Rep_ID) + (1|ExpID) + (1|Genotype:ExpID),
                 data = selection,
                 control=lmerControl(optimizer ="bobyqa",
                                     check.conv.singular = .makeCC(action = "message",  tol = 1e-4))) # define model
      random_effects <- ranef(model)
      BLUPs = random_effects$Genotype
      Genotypes = rownames(BLUPs)
      
      var.comp=VarCorr(model)
      Inter = fixef(model)
      Gvar =var.comp$Genotype 
      Expvar= var.comp$ExpID 
      GxEvar= var.comp$`Genotype:ExpID` 
      Repvar= var.comp$`Rep_ID:Genotype` 
      evar= attr(var.comp,'sc')^2
      H2= Gvar/(Gvar+Expvar+GxEvar+Repvar+evar)
    }
    
    cur_frame = data.frame(Genotype=Genotypes,
                           Trait=t_d_comb_HS$Trait[i],
                           DAT=t_d_comb_HS$DAT[i],
                           Value=BLUPs$`(Intercept)`)
    
    Var_comp_frame_HS[i,] = c(
      t_d_comb_HS$DAT[i],
      t_d_comb_HS$Trait[i],
      Inter,
      Gvar,
      Expvar,
      GxEvar,
      Repvar,
      evar,
      H2)
    if(i==1){
      BLUP_frame_HS = cur_frame
    }else{
      BLUP_frame_HS = rbind(BLUP_frame_HS,cur_frame)
    }
    
  }
  
  Var_comp_frame_HS <- Var_comp_frame_HS %>%
    mutate(WL=as.numeric(sub("VNIR_Plant_Spectrum_X","",Trait)))
  
  Var_comp_frame_HS$DAT = as.numeric(Var_comp_frame_HS$DAT)
  Var_comp_frame_HS$Intercept = as.numeric(Var_comp_frame_HS$Intercept)
  Var_comp_frame_HS$Gvar = as.numeric(Var_comp_frame_HS$Gvar)
  Var_comp_frame_HS$Expvar = as.numeric(Var_comp_frame_HS$Expvar)
  Var_comp_frame_HS$GxEvar = as.numeric(Var_comp_frame_HS$GxEvar)
  Var_comp_frame_HS$Repvar = as.numeric(Var_comp_frame_HS$Repvar)
  Var_comp_frame_HS$evar = as.numeric(Var_comp_frame_HS$evar)
  Var_comp_frame_HS$H2 = as.numeric(Var_comp_frame_HS$H2)
  
  
  
  write.csv(Var_comp_frame_HS,HSR_BLUP_variance_components_path,row.names = F)
  write.csv(BLUP_frame_HS,HSR_BLUP_path,row.names = F)

}else{
  Var_comp_frame_HS = read.csv(HSR_BLUP_variance_components_path) 
  BLUP_frame_HS = read.csv(HSR_BLUP_path) 
}

print(paste0("For ",sum(Var_comp_frame_HS$Gvar==0)," out of 2880 measurement instances, no genetic variance could be estimated."))

missing_set_HS = Var_comp_frame_HS[Var_comp_frame_HS$Gvar==0,]
table(missing_set_HS$DAT)
hist(x=missing_set_HS$DAT,,xlab = "DAT [d]")
hist(x=missing_set_HS$WL,xlab = "Wavelength [nm]")

## ---- Box-Cox transformation ----
## regular BLUPs
if(!(file.exists(BLUP_normalized_path)&file.exists(box_cox_parameter_file))){
parameter_frame = data.frame(DAT=NA,Trait=NA,mean=NA,sd=NA,lambda=NA)
BLUPs_normalized = data.frame(Genotype=NA,Trait=NA,DAT=NA,Value=NA)

for(i in 1:nrow(t_d_comb)){
  blup_selection <- BLUP_frame %>%
    filter(DAT == t_d_comb$DAT[i] & Trait == t_d_comb$Trait[i])
  if(nrow(blup_selection)==0){next}
  
  blups <- blup_selection$Value
  shifted <- blups + abs(min(blups)) + 1
  bc_obj = boxcox(shifted,standardize=T)
  
  parameter_frame[i,]=c(t_d_comb$DAT[i],t_d_comb$Trait[i],bc_obj$mean,bc_obj$sd,bc_obj$lambda)
  cur_frame = data.frame(Genotype=blup_selection$Genotype,
                         Trait=t_d_comb$Trait[i],
                         DAT=t_d_comb$DAT[i],
                         Value=bc_obj$x.t)
  if(i==1){
    BLUPs_normalized = cur_frame
  }else{
    BLUPs_normalized = rbind(BLUPs_normalized,cur_frame)
  }
}

write.csv(BLUPs_normalized,BLUP_normalized_path,row.names = F)
write.csv(parameter_frame,box_cox_parameter_file,row.names = F)

}else{
  BLUPs_normalized = read.csv(BLUP_normalized_path)
  parameter_frame = read.csv(box_cox_parameter_file)
}
## HSR BLUPs
if(!(file.exists(HSR_BLUP_normalized_path)&file.exists(HSR_box_cox_parameter_file))){
  HSR_parameter_frame = data.frame(DAT=NA,Trait=NA,mean=NA,sd=NA,lambda=NA)
  HSR_BLUPs_normalized = data.frame(Genotype=NA,Trait=NA,DAT=NA,Value=NA)
  for(i in 1:nrow(t_d_comb_HS)){
    blup_selection <- BLUP_frame_HS %>%
      filter(DAT == t_d_comb_HS$DAT[i] & Trait == t_d_comb_HS$Trait[i])
    if(nrow(blup_selection)==0){next}
    
    blups <- blup_selection$Value
    shifted <- blups + abs(min(blups)) + 1
    bc_obj = boxcox(shifted,standardize=T)
    
    HSR_parameter_frame[i,]=c(t_d_comb_HS$DAT[i],t_d_comb_HS$Trait[i],bc_obj$mean,bc_obj$sd,bc_obj$lambda)
    cur_frame = data.frame(Genotype=blup_selection$Genotype,
                           Trait=t_d_comb_HS$Trait[i],
                           DAT=t_d_comb_HS$DAT[i],
                           Value=bc_obj$x.t)
    if(i==1){
      HSR_BLUPs_normalized = cur_frame
    }else{
      HSR_BLUPs_normalized = rbind(HSR_BLUPs_normalized,cur_frame)
    }
  }
  
  write.csv(HSR_BLUPs_normalized,HSR_BLUP_normalized_path,row.names = F)
  write.csv(HSR_parameter_frame,HSR_box_cox_parameter_file,row.names = F)
  
}else{
  HSR_BLUPs_normalized = read.csv(HSR_BLUP_normalized_path)
  HSR_parameter_frame = read.csv(HSR_box_cox_parameter_file)
}

## ---- Heritability analysis general ----

# Average H2 & standard deviation
mean(Var_comp_frame$H2)
sd(Var_comp_frame$H2)

# Test heteroscedasticity of heritability btw. groups 
Var_comp_frame$Group = trait_groups$Group[match(Var_comp_frame$Trait,trait_groups$Trait)]
t=summary(aov(H2~Group,data=Var_comp_frame))


# Compare relative variance proportions (filter out single environment and single replicate)
Var_comp_frame_cut = Var_comp_frame[,c(4:8,10)] %>%
  filter(Expvar!=0 & Repvar!=0)
for(i in 1:5){
  component = colnames(Var_comp_frame_cut)[i]
  tot_var=rowSums(Var_comp_frame_cut[,1:5])
  rel_var=Var_comp_frame_cut[,i] / tot_var
  t=summary(aov(rel_var~Var_comp_frame_cut$Group))
  cur_row = data.frame(Component = component,
                       Mean_relative_variance = mean(rel_var,na.rm=T),
                       sd_relative_variance = sd(rel_var,na.rm=T),
                       F_value = t[[1]]$`F value`[1],
                       p_value = t[[1]]$`Pr(>F)`[1])
  if(i==1){
    relative_variance_frame = cur_row
  }else{
    relative_variance_frame = rbind(relative_variance_frame,cur_row)
  }
}
# significant p-values indicate heteroscedasticity between trait groups.
relative_variance_frame





## ---- Heritability analysis specific ----
## the following 6 data frames serve for general data exploration and summarize
## Group / Trait specific measures:
## mean | standard deviation | Coefficient of variation | Daily trend | adjusted R^2
## 
## The last two measures result from a linear fit of the data to the time series. 
## The adjusted R^2 value informs, how well a linear relationship with time 
## following the daily trend explains the observed data.


phenotype_trait_summary_statistics <- pheno_long %>%
  filter(!is.na(Value)) %>%
  group_by(Trait) %>%
  summarize(mean = mean(Value),
            sd = sd(Value),
            CoV = sd(Value)/mean(Value),
            dailyTrend = linearTemporalTrend2(DAT,Value)[[2]],
            Adj_R2 = linearTemporalTrend2(DAT,Value)[[3]]) %>%
  mutate(Group= trait_groups$Group[match(Trait,trait_groups$Trait)])

phenotype_group_summary_statistics <- pheno_long %>%
  filter(!is.na(Value)) %>%
  group_by(Group) %>%
  summarize(mean = mean(Value),
            sd = sd(Value),
            CoV = sd(Value)/mean(Value),
            dailyTrend = linearTemporalTrend2(DAT,Value)[[2]],
            Adj_R2 = linearTemporalTrend2(DAT,Value)[[3]])



BLUP_trait_summary_statistics <- BLUPs_normalized %>%
  filter(!is.na(Value)) %>%
  group_by(Trait) %>%
  summarize(mean = mean(Value),
            sd = sd(Value),
            CoV = sd(Value)/mean(Value),
            dailyTrend = linearTemporalTrend2(DAT,Value)[[2]],
            Adj_R2 = linearTemporalTrend2(DAT,Value)[[3]]) %>%
  mutate(Group= trait_groups$Group[match(Trait,trait_groups$Trait)])

BLUP_group_summary_statistics <- BLUPs_normalized %>%
  filter(!is.na(Value)) %>%
  group_by(Group) %>%
  summarize(mean = mean(Value),
            sd = sd(Value),
            CoV = sd(Value)/mean(Value),
            dailyTrend = linearTemporalTrend2(DAT,Value)[[2]],
            Adj_R2 = linearTemporalTrend2(DAT,Value)[[3]])


H2_trait_summary_statistics <- Var_comp_frame %>%
  mutate(Value=H2) %>%
  filter(!is.na(Value)) %>%
  group_by(Trait) %>%
  summarize(mean = mean(Value),
            sd = sd(Value),
            CoV = sd(Value)/mean(Value),
            dailyTrend = linearTemporalTrend2(DAT,Value)[[2]],
            Adj_R2 = linearTemporalTrend2(DAT,Value)[[3]]) %>%
  mutate(Group= trait_groups$Group[match(Trait,trait_groups$Trait)])

H2_group_summary_statistics <- Var_comp_frame %>%
  mutate(Value=H2) %>%
  filter(!is.na(Value)) %>%
  group_by(Group) %>%
  summarize(mean = mean(Value),
            sd = sd(Value),
            CoV = sd(Value)/mean(Value),
            dailyTrend = linearTemporalTrend2(DAT,Value)[[2]],
            Adj_R2 = linearTemporalTrend2(DAT,Value)[[3]])



## Correlation of heritability estimates for traits in each group across time points.
# i.e. do the heritability estimates of two traits in the same group show similarities 
# in value distribution across time.
nf=T
for(group in unique(trait_groups$Group)){
  H2_cor_mat = Var_comp_frame %>%
    filter(Group==group) %>%
    dplyr::select(Trait,DAT,H2) %>%
    pivot_wider(id_cols = DAT,names_from = Trait,values_from = H2)
  
  c=cor(x=as.matrix(H2_cor_mat[,-1]))
  current = data.frame(Group = group,
                       corr_mean =  mean(c[upper.tri(c)],na.rm=T),
                       corr_sd = sd(c[upper.tri(c)],na.rm=T))
  if(nf){
    corr_frame = current
    nf=F
  }else{
    corr_frame = rbind(PS_corr_frame,current)
  }
}





## ---- Photosynthesis traits ----
Photosynthesis_TraitGroups = c("CC","CL","HL","HL.LL")
Photosynthesis_Heritability = Var_comp_frame %>% 
  filter(Group %in% Photosynthesis_TraitGroups)

## Average heritability in photosynthesis trait groups

summary=rbind(data.frame(Group="All PS",
                         H2_mean=mean(Photosynthesis_Heritability$H2),
                         H2_sd=sd(Photosynthesis_Heritability$H2)))
summary

# Heritability across time
H2_group_summary_statistics %>%
  filter(Group %in% Photosynthesis_TraitGroups)

# Correlation among heritability estimates through time
corr_frame %>%
  filter(Group %in% Photosynthesis_TraitGroups)

## Latent factor analysis
# Dimensionality reduction (PCA, with a cumulative variance threshold of 90%) on
# phenotype (normalized BLUPs) to find fewer factors explaining large 
# proportions of variance
BLUPs_normalized$Group = trait_groups$Group[match(BLUPs_normalized$Trait,trait_groups$Trait)]

lat_fac_df <- BLUPs_normalized %>%
  filter(Group %in% c("CC","CL","HL","HL.LL")) %>%
  mutate(ID=paste0(Genotype,"_",DAT)) %>%
  dplyr::select(ID,Group,Trait,Value)

lat_fac_df_CC <- lat_fac_df %>%
  filter(Group =="CC") %>%
  pivot_wider(id_cols = ID,names_from = Trait,values_from = Value)

lat_fac_df_CL <- lat_fac_df %>%
  filter(Group =="CL") %>%
  pivot_wider(id_cols = ID,names_from = Trait,values_from = Value)

lat_fac_df_HL <- lat_fac_df %>%
  filter(Group =="HL") %>%
  pivot_wider(id_cols = ID,names_from = Trait,values_from = Value)

lat_fac_df_HLLL <- lat_fac_df %>%
  filter(Group =="HL.LL") %>%
  pivot_wider(id_cols = ID,names_from = Trait,values_from = Value)

pca_reduction(lat_fac_df_CC[,-1],0.9)
pca_reduction(lat_fac_df_CL[,-1],0.9)
pca_reduction(lat_fac_df_HL[,-1],0.9)
pca_reduction(lat_fac_df_HLLL[,-1],0.9)



## ---- VNIR traits ----
H2_group_summary_statistics %>%
  filter(Group=="VNIR")

## Correlation of heritability estimates for traits in each group across time points.
corr_frame[corr_frame$Group=="VNIR",]



## ---- IR traits ----
## Correlation between the 2 traits across time.
corr_frame[corr_frame$Group=="IR",]


## Comparison average plant temperature / delta temperature
IR_plot = Var_comp_frame %>%
  filter(Group == "IR") %>%
  dplyr::select(Trait,DAT,H2)

ggplot(IR_plot,aes(x=DAT,y=H2,color=Trait))+
  geom_line() +
  theme_linedraw()

IR_comp_frame = IR_plot %>%
  group_by(DAT) %>%
  summarize(Delta = H2[Trait=="IR1.Probes_Plant_deltaT"] - H2[Trait=="IR1_Plant_Temp.avg"])

print(paste0("Average H2 difference between delta temperature and average plant temperature is: ",
             round(mean(IR_comp_frame$Delta),2),
             "±",
             round(sd(IR_comp_frame$Delta),2)))

      
CV_pat = pheno_long %>%
  filter(Trait=="IR1_Plant_Temp.avg") %>%
  group_by(DAT) %>%
  summarize(CoV = sd(Value) / mean(Value))

print(paste0("Average coefficient of variation of plant average temperature is:",
             round(mean(CV_pat$CoV),2),
             "±",
             round(sd(CV_pat$CoV),2)))

pheno_dT = pheno_long %>%
  filter(Trait=="IR1.Probes_Plant_deltaT")%>%
  group_by(DAT) %>%
  summarise(mean=mean(Value),sd=sd(Value))

t = linearTemporalTrend2(pheno_dT$DAT,pheno_dT$mean)
print(paste0("temperature regulation increased ",round(t[[2]],3)," °C per day"))

## ---- RGB traits ----
H2_group_summary_statistics %>%
  filter(Group=="RGB")

H2_trait_summary_statistics %>% 
  filter(Group=="RGB")

## Correlation of heritability estimates for traits in each group across time points.
corr_frame[corr_frame$Group=="RGB",]


## ---- SW traits ----
H2_group_summary_statistics %>%
  filter(Group=="SW")

H2_trait_summary_statistics %>% 
  filter(Group=="SW")

## Correlation of heritability estimates for traits in each group across time points.
corr_frame[corr_frame$Group=="SW",]

## Coefficient of variation in measured phenotype
phenotype_trait_summary_statistics %>%
  filter(Group == "SW")

## ---- Harvest traits ---- ####
H2_group_summary_statistics %>%
  filter(Group=="H")

H2_trait_summary_statistics %>% 
  filter(Group=="H")




## ---- Set of measurement insances with zero genotypic variance ----
missing_set = BLUPs_normalized %>%
  filter(is.na(Value)) %>%
  dplyr::select(Trait,DAT,Group) %>%
  distinct(.keep_all=T)


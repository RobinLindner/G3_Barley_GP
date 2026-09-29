## ---------------------------
##
## Script name: 4_GenomicPrediction
##
## Purpose of script: GenomicPrediction with the four discussed methods: GBLUP,
##                    HBLUP, G+HBLUP, and MegaLMM
##
## Author: M.Sc. Robin Lindner
##
## Date Created: 2025-07-10
##
## Copyright (c) Robin Lindner, 2025
## Email: robin.lindner@uni-potsdam.de
##
## ---------------------------
##
## Notes: 
## - Univariate GP requires large amounts of RAM. Can be outsourced.
## - MegaLMM GP requires large amounts of computation time and storage. Should be outsourced & Parallelized.
##
## ---------------------------

## set environment

source("0_utils.R") 

## ---------------------------

## ---- create a directory for GP results ---- 
dir.create("../GP_Results")
outpath = "../GP_Results"

## ---- load GP input data ----

focal_traits = c("RGB1_Plant_Avg_HEIGHT_MM","VNIR_Plant_NDVI.avg")

# Read Kinship matrix as covariance strukture K
K = read.csv(GRM_path,row.names = 1)

# Genotypes that had measurements for all time points & traits & HSR data
genotypes = read.table(geno4GP_file)$X

# Read numeric genotype to extract marker fixed effects (MFE)
marker_geno = read.table(numeric_geno_file)
marker_sites = read.table(sites_file,header=T)
marker_taxa = read.table(taxa_file,header=T)
dimnames(marker_geno)=list(x=marker_sites$Name,y=marker_taxa$Taxa)

# mxn => nxm
marker_geno= t(marker_geno)

# Use significant associations to select MFE
sig_associations = read.csv(sig_associations_file)
assoc_cut <- sig_associations %>%
  filter(Trait %in% focal_traits) %>%
  select(SNP,Value,DAT,Trait) %>%
  distinct(.keep_all=T)

# Prepare the CV matrix for reproducibility
nfold=5
nrun=20
CV_mat = replicate(nrun,partition(1:length(genotypes),nfold))
rownames(CV_mat) = genotypes
colnames(CV_mat) = paste0("run_",c(1:nrun))

# write out the CV matrix for use in MegaLMM.
write.csv(CV_mat,GP_CV_matrix_file)


# Prepare the normalized BLUPs (i.e. focal traits)
BLUPs_normalized = read.csv(BLUP_normalized_path) %>%
  filter(Genotype %in% genotypes)
names(BLUPs_normalized)[c(1,4)]=c("X","BLUP")

# Prepare the normalized BLUPs of the HSR data (i.e. secondary traits)
HSR_BLUPs_normalized = read.csv(HSR_BLUP_normalized_path)


## ---- Univariate GP ----

all_trait_acc_df=data.frame()
all_trait_pred_df=data.frame()
all_trait_MFE_acc_df=data.frame()
all_trait_MFE_pred_df=data.frame()
all_trait_MFE_fe_df=data.frame()

dats_HS = c(14,21,28,35,42)

for(i in 1:length(focal_traits)){
  
  trait = focal_traits[i]
  if(i==1){
    dats_foc = dats_HS + 1
  }else{
    dats_foc = dats_HS
  }
  
  trait_acc_df=data.frame()
  trait_pred_df=data.frame()
  for(j in 1:length(dats_foc)){
    dat_foc=dats_foc[j]
    
    sig_snp <- assoc_cut %>%
      filter(Trait==trait,DAT==dat_foc) %>%
      select(SNP)
    
    X <- as.data.frame(marker_geno[,sig_snp$SNP])
    names(X)=sig_snp$SNP
    
    dat_HS_blups <- all_BLUP_HS %>%
      filter(DAT==dats_HS[j])
    
    HS_mat = pivot_wider(dat_HS_blups,id_cols = c(X),names_from = Trait,values_from = BLUP)%>%
      select(where(~any(. !=1))) %>%
      column_to_rownames(var="X")
    
    # prevent singularity
    rem=findLinearCombos(as.matrix(HS_mat))$remove
    HS_mat=HS_mat[-rem]
    
    H=as.matrix(HS_mat) %*% t(as.matrix(HS_mat)) / ncol(HS_mat) 
    
    print(paste0("DAT:",j," Trait:",i))
    print(paste0("number of sigificant SNPs: ",ncol(X)))
    
    # Returns a list containing 
    # $Accuracy:    | GBLUP Accuracy | HBLUP Accuracy | G+HBLUP Accuracy | Run | Fold
    # $Predictions: | Genotype | BLUPs | GBLUP Prediction | HBLUP Prediction | G+HBLUP Prediction | Test/Train | Run | Fold |
    res = nFoldCV_lm_combined(all_BLUPs,trait, dat_foc,K, H, CV_mat ,genotypes)
    
    # Returns a list containing 
    # $Accuracy:    | GBLUP Accuracy | HBLUP Accuracy | G+HBLUP Accuracy | GBLUP Accuracy NonAdj | HBLUP Accuracy NonAdj | G+HBLUP Accuracy NonAdj | Run | Fold | #Fixed effects | #reduced columns in genotype|
    # $Predictions: | Genotype | BLUPs | GBLUP Prediction | HBLUP Prediction | G+HBLUP Prediction | Test/Train | Run | Fold |
    # $FE_sizes:    | Fixed effect ID | GBLUP effect size | HBLUP effect size | G+HBLUP effect size | Run | Fold |
    # NonAdj: prediction accuracy without considering marker effects cor(u,y) <=> cor(y',y)
    res_MFE = nFoldCV_lm_combined_MFE(all_BLUPs,trait, dat_foc, X ,K, H, CV_mat, genotypes)
    
    if(i==1 & j==1){
      trait_acc_df = cbind(res$Accuracy,data.frame(DAT=rep(dat_foc,nrow(res$Accuracy))))
      trait_pred_df = cbind(res$Predictions,data.frame(DAT=rep(dat_foc,nrow(res$Predictions))))
      trait_MFE_acc_df = cbind(res_MFE$Accuracy,data.frame(DAT=rep(dat_foc,nrow(res_MFE$Accuracy))))
      trait_MFE_pred_df = cbind(res_MFE$Predictions,data.frame(DAT=rep(dat_foc,nrow(res_MFE$Predictions))))
      trait_MFE_fe_df = cbind(res_MFE$FE_sizes,data.frame(DAT=rep(dat_foc,nrow(res_MFE$FE_sizes))))
    }else{
      trait_acc_df = rbind(trait_acc_df,cbind(res$Accuracy,data.frame(DAT=rep(dat_foc,nrow(res$Accuracy)))))
      trait_pred_df = rbind(trait_pred_df,cbind(res$Predictions,data.frame(DAT=rep(dat_foc,nrow(res$Predictions)))))
      trait_MFE_acc_df = rbind(trait_acc_df,cbind(res_MFE$Accuracy,data.frame(DAT=rep(dat_foc,nrow(res_MFE$Accuracy)))))
      trait_MFE_pred_df = rbind(trait_pred_df,cbind(res_MFE$Predictions,data.frame(DAT=rep(dat_foc,nrow(res_MFE$Predictions)))))
      trait_MFE_fe_df = rbind(trait_fe_df,cbind(res_MFE$FE_sizes,data.frame(DAT=rep(dat_foc,nrow(res_MFE$FE_sizes)))))
      
    }
  }
  write.csv(trait_acc_df,paste0(out_path,"/",trait,"_accuracy.csv"),row.names = F)
  write.csv(trait_pred_df,paste0(out_path,"/",trait,"_predictions.csv"),row.names = F)
  write.csv(trait_MFE_acc_df,paste0(out_path,"/",trait,"_MFE_accuracy.csv"),row.names = F)
  write.csv(trait_MFE_pred_df,paste0(out_path,"/",trait,"_MFE_predictions.csv"),row.names = F)
  write.csv(trait_MFE_fe_df,paste0(out_path,"/",trait,"_MFE_fixed_effects.csv"),row.names = F)
  
  if(i==1){
    all_trait_acc_df = cbind(trait_acc_df,data.frame(Trait=rep(trait,nrow(trait_acc_df))))
    all_trait_pred_df = cbind(trait_pred_df,data.frame(Trait=rep(trait,nrow(trait_pred_df))))
    all_trait_MFE_acc_df = cbind(trait_MFE_acc_df,data.frame(Trait=rep(trait,nrow(all_trait_MFE_acc_df))))
    all_trait_MFE_pred_df = cbind(trait_MFE_pred_df,data.frame(Trait=rep(trait,nrow(all_trait_MFE_pred_df))))
    all_trait_MFE_fe_df = cbind(trait_MFE_fe_df,data.frame(Trait=rep(trait,nrow(all_trait_MFE_fe_df))))
    
  }else{
    all_trait_acc_df = rbind(all_trait_acc_df,cbind(trait_acc_df,data.frame(Trait=rep(trait,nrow(trait_acc_df)))))
    all_trait_pred_df = rbind(all_trait_pred_df,cbind(trait_pred_df,data.frame(Trait=rep(trait,nrow(trait_pred_df)))))
    all_trait_MFE_acc_df = rbind(all_trait_MFE_acc_df,cbind(trait_MFE_acc_df,data.frame(Trait=rep(trait,nrow(all_trait_MFE_acc_df)))))
    all_trait_MFE_pred_df = rbind(all_trait_MFE_pred_df,cbind(trait_MFE_pred_df,data.frame(Trait=rep(trait,nrow(all_trait_MFE_pred_df)))))
    all_trait_MFE_fe_df = rbind(all_trait_MFE_fe_df,cbind(trait_MFE_fe_df,data.frame(Trait=rep(trait,nrow(all_trait_MFE_fe_df)))))
    
  }
}

write.csv(all_trait_acc_df,paste0(out_path,"/all_trait_SV_noMFE_accuracy.csv"))
write.csv(all_trait_pred_df,paste0(out_path,"/all_trait_SV_noMFE_predictions.csv"))
write.csv(all_trait_MFE_acc_df,paste0(out_path,"/all_trait_SV_MFE_accuracy.csv"))
write.csv(all_trait_MFE_pred_df,paste0(out_path,"/all_trait_SV_MFE_predictions.csv"))
write.csv(all_trait_MFE_fe_df,paste0(out_path,"/all_trait_SV_MFE_fixed_effects.csv"))

## ---- MegaLMM ----

## For this computation a HPC is recommended, since the storage demand and computing 
## time are very high.

# Computations can be run in parallel by
# 
# for (trait) {
# for (time point) {
# for (run) {
# 
# Rscript MegaLMM_GP.R [BLUPs_normalized_path] [HSR_BLUPs_normalized_path] [GRM_path] [GPgenotypes_file] [GP_CV_matrix_file] [output_folder] [trait] [dat] [run] 
# 
# Rscript MegaLMM_MFE_GP.R [BLUPs_normalized_path] [HSR_BLUPs_normalized_path] [GRM_path] [GPgenotypes_file] [sig_associations_file] [numeric_geno_file] [GP_CV_matrix_file] [output_folder] [trait] [dat] [run] 
# }
# }
# }


# ---- Loading data ----
files = list.files("../Data/Generated/GenomicPrediction/UVGP/")

UV_accuracy_df = data.frame()
UV_prediction_df = data.frame()
UV_FE_df = data.frame()

acc_nf = T
pred_nf= T
fe_nf = T
for(file in files){
  if(grepl("_accuracy.csv",file)){
    if(grepl("_MFE_",file)){
      scenario = sub("_MFE_accuracy.csv","",file)
    }else{
      scenario = sub("_accuracy.csv","",file)
    }
    t = strsplit(scenario,"_")
    dat = t[[1]][length(t[[1]])]
    trait = paste0((t[[1]][-length(t[[1]])]),collapse = "_")
    df = read.csv(paste0("../Data/Generated/GenomicPrediction/UVGP_old//",file))
    df$Trait = trait
    df$Dat = dat
    df$MFE = grepl("_MFE_",file)
    if(acc_nf){
      UV_accuracy_df = df
      acc_nf = F
    }else{
      UV_accuracy_df = bind_rows(UV_accuracy_df, df)
    }
  }
  else if(grepl("_predictions.csv", file)){
    if(grepl("_MFE_", file)){
      scenario = sub("_MFE_predictions.csv","",file)
    }else{
      scenario = sub("_predictions.csv","",file)
    }
    t = strsplit(scenario,"_")
    dat = t[[1]][length(t[[1]])]
    trait = paste0((t[[1]][-length(t[[1]])]),collapse = "_")
    df = read.csv(paste0("../Data/Generated/GenomicPrediction/UVGP/",file))
    df$Trait = trait
    df$Dat = dat
    df$MFE = grepl("_MFE_", file)
    if(pred_nf){
      UV_prediction_df = df
      pred_nf = F
    }else{
      UV_prediction_df = bind_rows(UV_prediction_df, df)
    }
  }
  else if(grepl("_fixed_effects", file)){
    scenario = sub("_MFE_fixed_effects.csv","",file)
    t = strsplit(scenario,"_")
    dat = t[[1]][length(t[[1]])]
    trait = paste0((t[[1]][-length(t[[1]])]),collapse = "_")
    df = read.csv(paste0("../Data/Generated/GenomicPrediction/UVGP/",file))
    df$Trait = trait
    df$Dat = dat
    if(fe_nf){
      UV_FE_df = df
      fe_nf = F
    }else{
      UV_FE_df = bind_rows(UV_FE_df, df)
    }
  }
}

acc_files = list.files("../Data/Generated/GenomicPrediction/MVGP/Accuracy/")
pred_files = list.files("../Data/Generated/GenomicPrediction/MVGP/Predictions/")

MV_accuracy_df = data.frame()
MV_prediction_df = data.frame()
MV_FE_df = data.frame()

acc_nf = T
pred_nf= T
fe_nf = T

for(file in acc_files){
  if(grepl("_CV2_",file)){
    CV2 = T
  }else{
    CV2 = F
  }
  if(grepl("_MFE_",file)){
    MFE = T
  }else{
    MFE = F
  }
  if(MFE){
    if(CV2){
      scenario = sub("MFE_CV2_accuracy.csv","",file)
    }else{
      scenario = sub("MFE_accuracy.csv","",file)
    }
    
  }else{
    if(CV2){
      scenario = sub("_CV2_accuracy.csv","",file)
    }else{
      scenario = sub("_accuracy.csv","",file)
    }
    
  }
  t = strsplit(scenario,"_")
  dat = t[[1]][length(t[[1]])]
  trait = paste0((t[[1]][-length(t[[1]])]),collapse = "_")
  df = read.csv(paste0("../Data/Generated/GenomicPrediction/MVGP/Accuracy/",file))
  #if(MFE){
  #  colnames(df)[c(3,2)] = c("Uhat_accuracy","Eta_mean_accuracy")
  #}else{
  #  colnames(df)[c(2,3)] = c("MegaLMM","MegaLMM_secondary")
  #}
  df$MFE = MFE
  df$CV2 = CV2
  if(acc_nf){
    MV_accuracy_df = df
    acc_nf = F
  }else{
    MV_accuracy_df = bind_rows(MV_accuracy_df, df)
  }
}
for(file in pred_files){
  if(grepl("_prediction", file)){
    if(grepl("_CV2_",file)){
      CV2 = T
    }else{
      CV2 = F
    }
    if(grepl("_MFE_",file)){
      MFE = T
    }else{
      MFE = F
    }
    if(MFE){
      if(CV2){
        scenario = sub("MFE_CV2_predicion.csv","",file)
      }else{
        scenario = sub("MFE_predicion.csv","",file)
      }
      
    }else{
      if(CV2){
        scenario = sub("_CV2_predicion.csv","",file)
      }else{
        scenario = sub("_predicion.csv","",file)
      }
      
    }
    t = strsplit(scenario,"_")
    dat = t[[1]][length(t[[1]])]
    trait = paste0((t[[1]][-length(t[[1]])]),collapse = "_")
    df = read.csv(paste0("../Data/Generated/GenomicPrediction/MVGP/Predictions/",file),row.names = 1)
    df$MFE = MFE
    df$CV2 = CV2
    if(pred_nf){
      MV_prediction_df = df
      pred_nf = F
    }else{
      MV_prediction_df = bind_rows(MV_prediction_df, df)
    }
  }
  else if(grepl("fixedEffects", file)){
    if(grepl("_CV2_",file)){
      CV2 = T
    }else{
      CV2 = F
    }
    if(CV2){
      scenario = sub("MFE_CV2_fixedEffects.csv","",file)
    }else{
      scenario = sub("MFE_fixedEffects.csv","",file)
    }
    t = strsplit(scenario,"_")
    dat = t[[1]][length(t[[1]])]
    trait = paste0((t[[1]][-length(t[[1]])]),collapse = "_")
    df = read.csv(paste0("../Data/Generated/GenomicPrediction/MVGP/Predictions/",file),row.names = 1)
    df$Trait = trait
    df$Dat = dat
    df$CV2 = CV2
    if(fe_nf){
      MV_FE_df = df
      fe_nf = F
    }else{
      MV_FE_df = bind_rows(MV_FE_df, df)
    }
  }
}

MV_FE_df$Model = NA
MV_FE_df$Model[MV_FE_df$CV2] = "MegaLMM"
MV_FE_df$Model[!MV_FE_df$CV2] = "MegaLMM CV1"

MV_FE_df_c = MV_FE_df %>%
  dplyr::select(SNPs,rep,Trait,DAT,Model,FE)

colnames(MV_FE_df_c) = c("EffectTerm","Repetition","Trait","DAT","Model","Value")

UV_FE_df_c=UV_FE_df %>% 
  pivot_longer(cols=c(GBLUP_fe,HBLUP_fe,GHBLUP_fe),values_to = "Value",names_to = "Model") %>%
  mutate(Model = sub("_fe","",Model)) %>%
  dplyr::select(FE,Repetition,Trait,Dat,Model,Value)

colnames(UV_FE_df_c) = c("EffectTerm","Repetition","Trait","DAT","Model","Value")

FE_df = rbind(UV_FE_df_c,MV_FE_df_c)
write.csv(FE_df,"../Data/Generated/GenomicPrediction/Fixed_effects_data.csv",row.names = F)


AccuracyDf = read.csv("../Figures/GP_accuracy_data.csv",row.names = 1)


# ---- Evaluation ----

summarize_trait_dat = AccuracyDf %>%
  group_by(Trait,Dat) %>%
  summarize(Mean=mean(CorrectedPA),
            sd = sd(CorrectedPA))


lm(data=AccuracyDf[AccuracyDf$Trait=="RGB1_Plant_Avg_HEIGHT_MM",],formula = CorrectedPA~Dat)
lm(data=AccuracyDf[AccuracyDf$Trait=="VNIR_Plant_NDVI.avg",],formula = CorrectedPA~Dat)

temp = merge(AccuracyDf,h2_df)
cor(temp$CorrectedPA,temp$h2,method = "spearman")

nf=T
models = unique(AccuracyDf$Model)
for(i in 1:4){
  for(j in (i+1):5){
    for(dat in unique(AccuracyDf$Dat)){
      for(trait in unique(AccuracyDf$Trait)){
        model1=models[i]
        model2=models[j]
         t1=AccuracyDf %>%
          filter(Model==model1 & Dat ==dat & !MFE & Trait == trait)
         t2=AccuracyDf %>%
           filter(Model==model2 & Dat ==dat & !MFE & Trait == trait)
         temp = t.test(t1$CorrectedPA,t2$CorrectedPA,)
         cur = data.frame(Trait=trait,
                          DAT=dat,
                          Model1 = model1,
                          Model2 = model2,
                          t=temp$statistic,
                          p_val = temp$p.value)
         if(nf){
           significance_testing_df = cur
           nf=F
         }else{
           significance_testing_df=rbind(significance_testing_df,cur)
         }
      }
    }
  }
}

summarize_model_trait = AccuracyDf %>%
  group_by(Model,Trait) %>%
  summarize(Mean=mean(CorrectedPA),
            sd = sd(CorrectedPA))

nf=T
models = unique(AccuracyDf$Model)
for(i in 1:5){
  for(dat in unique(AccuracyDf$Dat)){
    for(trait in unique(AccuracyDf$Trait)){
      model1=models[i]
      model2=models[j]
      t1=AccuracyDf %>%
        filter(Model==models[i] & Dat ==dat & !MFE & Trait == trait)
      t2=AccuracyDf %>%
        filter(Model==models[i] & Dat ==dat & MFE & Trait == trait)
      diff = t1$CorrectedPA[order(t1$Repetition)] - t2$CorrectedPA[order(t2$Repetition)]
      diff_rel = t2$CorrectedPA[order(t1$Repetition)] / t1$CorrectedPA[order(t2$Repetition)]
      temp = t.test(t1$CorrectedPA,t2$CorrectedPA)
      cur = data.frame(Trait=trait,
                       DAT=dat,
                       Model1 = models[i],
                       t=temp$statistic,
                       p_val = temp$p.value,
                       diff_mean = mean(diff),
                       diff_sd = sd(diff),
                       rel_mean = mean(diff_rel),
                       rel_sd = sd(diff_rel))
      if(nf){
        mfe_significance_testing_df = cur
        nf=F
      }else{
        mfe_significance_testing_df=rbind(mfe_significance_testing_df,cur)
      }
    }
    
  }
}

mfe_significance_testing_df$sig=mfe_significance_testing_df$p_val<0.05

nf=T
models = unique(AccuracyDf$Model)
traits = unique(AccuracyDf$Trait)
for(i in 1:5){
  for(dat in unique(AccuracyDf$Dat)){
    for(MFE in c(T,F)){
      model1=models[i]
      model2=models[j]
      t1=AccuracyDf %>%
        filter(Model==models[i] & Dat ==dat & !MFE & Trait == traits[1])
      t2=AccuracyDf %>%
        filter(Model==models[i] & Dat ==dat & MFE & Trait == traits[2])
      diff = t1$CorrectedPA[order(t1$Repetition)] - t2$CorrectedPA[order(t2$Repetition)]
      
      temp = t.test(t1$CorrectedPA,t2$CorrectedPA,)
      cur = data.frame(MFE=MFE,
                       DAT=dat,
                       Model1 = models[i],
                       t=temp$statistic,
                       p_val = temp$p.value,
                       diff_mean = mean(diff),
                       diff_sd = sd(diff))
      if(nf){
        trait_significance_testing_df = cur
        nf=F
      }else{
        trait_significance_testing_df=rbind(trait_significance_testing_df,cur)
      }
    }
    
  }
}

tt = trait_significance_testing_df %>%
  filter(!MFE)


mfe_significance_testing_df$sig=mfe_significance_testing_df$p_val<0.05





ph_acc = temp %>%
  filter(Trait == "RGB1_Plant_Avg_HEIGHT_MM")
ph_acc = ph_acc[order(ph_acc$Dat,ph_acc$Repetition,ph_acc$Model),]

ph_acc_mfe = ph_acc%>%
  filter(MFE)
ph_acc_no_mfe = ph_acc%>%
  filter(!MFE)
diff_vec = ph_acc_mfe$CorrectedPA - ph_acc_no_mfe$CorrectedPA

mean(diff_vec)
sd(diff_vec)

cor(ph_acc_mfe$h2,diff_vec)

temp = mfe_significance_testing_df[mfe_significance_testing_df$Model1!="MegaLMM_CV1",]


write.csv(mfe_significance_testing_df[,c(1,2,3,4,5)],"../Supplements/MFE_significance_testing.csv",row.names = F)

snp_counts = FE_df %>% 
  group_by(Trait,DAT) %>%
  summarize(nSNP = length(unique(EffectTerm))-1)

FE_df_wide = FE_df %>%
  pivot_wider(values_from = Value,names_from = Model)

cor_mat = cor(FE_df_wide[5:9],use = "p")
mean(cor_mat[1:3,4])
sd(cor_mat[1:3,4])

CoV <- function(x){sd(x)/mean(x)}

cov_stats=FE_df_wide %>% 
  group_by(EffectTerm,Trait,DAT) %>%
  summarize(GBLUP = abs(CoV(GBLUP)),
            HBLUP = abs(CoV(HBLUP)),
            GHBLUP = abs(CoV(GHBLUP)),
            MegaLMM = abs(CoV(MegaLMM)),
            `MegaLMM CV1` = abs(CoV(`MegaLMM CV1`))) %>%
  filter(Trait!="SC_Plant_Weight")



FE_strongest_effects_table = FE_df_wide %>%
  filter(Trait!="SC_Plant_Weight") %>%
  filter(EffectTerm!="Intercept") %>%
  group_by(EffectTerm,Trait,DAT) %>%
  summarize(GBLUP = mean(GBLUP),
            HBLUP = mean(HBLUP),
            GHBLUP = mean(GHBLUP),
            MegaLMM = mean(MegaLMM),
            `MegaLMM CV1` = mean(`MegaLMM CV1`))

write.csv(FE_strongest_effects_table,"../Supplements/FE_average_effect_sizes_table.csv",row.names = F)


rel_diff_frame =AccuracyDf %>%
  dplyr::select(MFE,Dat,Trait,Model,CorrectedPA,Repetition) %>%
  pivot_wider(values_from = CorrectedPA,names_from = MFE)
colnames(rel_diff_frame)[c(5,6)] = c("noMFE","MFE")
rel_diff_frame = rel_diff_frame%>%
  group_by(Trait,Model) %>%
  summarize(MFE = mean(MFE),
            noMFE = mean(noMFE)) %>%
  mutate(RelativeDifference = (MFE-noMFE)/abs(noMFE))

ph_rd_frame = rel_diff_frame %>%
  group_by(Trait,Model,Dat) %>%
  summarize(mean(RelativeDifference),
            sd(RelativeDifference))

mean(rel_diff_frame$RelativeDifference[rel_diff_frame$Trait==traits[1]])

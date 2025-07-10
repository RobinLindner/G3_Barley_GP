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

focal_traits = c("RGB1_Plant_Avg_HEIGHT_MM","VNIR_Plant_NDVI.avg","SC_Plant_Weight")

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
for(i in 1:length(focal_traits)){
  
  trait = focal_traits[i]
  dats_HS = c(14,21,28,35,42)
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





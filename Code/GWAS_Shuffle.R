## ---------------------------
##
## Script name: GWAS_Shuffle.R
##
## Purpose of script: Find training/testing partitions that share SNP for each trait & dat
##
## Author: M.Sc. Robin Lindner
##
## Date Created: 2026-02-03
##
## Copyright (c) Robin Lindner, 2026
## Email: robin.lindner@uni-potsdam.de
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac

setwd("~/Documents/Arbeit/")    

## ---------------------------

## load up the packages we will need:  
source("0_utils.R")
## ---------------------------

focal_traits = c("RGB1_Plant_Avg_HEIGHT_MM","VNIR_Plant_NDVI.avg","SC_Plant_Weight")
dats_HS = c(14,21,28,35,42)
K = read.csv(GRM_path,row.names = 1)
gt = read.table(geno4GP_file)$X
marker_geno = read.table(numeric_geno_file)
marker_sites = read.table(sites_file,header=T)
marker_taxa = read.table(taxa_file,header=T)
dimnames(marker_geno)=list(x=marker_sites$Name,y=marker_taxa$Taxa)


nfold=10
nrun=100
CV_mat = replicate(nrun,partition(1:length(gt),nfold))
rownames(CV_mat) = gt
colnames(CV_mat) = paste0("run_",c(1:nrun))

write.csv(CV_mat,"../Supplements/GP_CV_mat2.csv")
# mxn => nxm
marker_geno= t(marker_geno)

full_df = data.frame()
nf=T
for(t in 1:length(focal_traits)){
  trait = focal_traits[t]
  for(d in 1:length(dats_HS)){
    dat = dats_HS[d]
    BLUPS_foc_spec <- all_BLUPs %>%
      filter(Trait == trait) %>%
      filter(DAT == dat) %>%
      filter(X %in% gt) %>%
      dplyr::select(X,BLUP)
    for(i in 1:n_runs){
      fold_vec=CV_mat[,i]
      for(j in 1:n_folds){
        test_idx=which(fold_vec==j)
        train_idx = which(fold_vec!=j)
        
        geno_train = marker_geno[BLUPS_foc_spec$X[train_idx],]
        K_train = K[BLUPS_foc_spec$X[train_idx],BLUPS_foc_spec$X[train_idx]]
        Y_train = BLUPS_foc_spec[train_idx,]
        colnames(Y_train) = c("Taxa",trait)
        sig_snp = GetSignificantAssociationsForTrainingSet(Y_train = Y_train,
                                                           geno_train =  geno_train,
                                                           K_train = K_train)
        if(length(sig_snp)==0){
          sig_snp = NA
        }
        cur_df = data.frame(SNP_idxs = sig_snp,
                            Trait= trait,
                            DAT = DAT,
                            run = i,
                            fold = j)
        
        if(nf){
          full_df = cur_df
          nf=F
        }else{
          full_df = rbind(full_df,cur_df)
        }
      }
    }
  }
}


    





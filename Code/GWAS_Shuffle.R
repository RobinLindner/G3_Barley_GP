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
## Notes: Run on the cluster via 
##   sbatch --array=0-749%20 Code/GWAS_shuffle.sh
##   sbatch --array=750-1499%20 Code/GWAS_shuffle.sh
## ---------------------------

## set working directory for Mac
setwd("Code")
## ---------------------------

## load up the packages we will need:  
## ---------------------------
source("0_utils.R")

args = commandArgs(trailingOnly = T)

# CV matrix
CV_mat_file = paste0("../",args[1])

# output folder
out_path = paste0("../", args[2])

# trait
trait = args[3]

# DAT
dat = as.numeric(args[4])

# runID
run = as.numeric(args[5])

## parallel access to files is staggered to prevent reading errors 
read_csv_retry <- function(file,rn=F, max_retries = 3, delay = 5) {
  attempts <- 0
  while (attempts < max_retries) {
    attempts <- attempts + 1
    tryCatch({
      if(rn){
        data <- read.csv(file,row.names = 1)
      }else{
        data <- read.csv(file)
      }
      return(data)
    }, error = function(e) {
      if (attempts == max_retries) stop("Max retries reached. Unable to read file.")
      Sys.sleep(delay) # Wait before retrying
    })
  }
}

read_table_retry <- function(file,header=F, max_retries = 3, delay = 5) {
  attempts <- 0
  while (attempts < max_retries) {
    attempts <- attempts + 1
    tryCatch({
      data <- read.table(file,header=header)
      return(data)
    }, error = function(e) {
      if (attempts == max_retries) stop("Max retries reached. Unable to read file.")
      Sys.sleep(delay) # Wait before retrying
    })
  }
}



K = read_csv_retry(GRM_path,rn = T)
print("K loaded")
gt = read_table_retry(geno4GP_file,header=T)$X
print("gt loaded")
marker_geno = read_table_retry(numeric_geno_file)
print("genotypes loaded")
marker_sites = read_table_retry(sites_file,header=T)
print("marker sites loaded")
marker_taxa = read_table_retry(taxa_file,header=T)
print("marker taxa loaded")
dimnames(marker_geno)=list(x=marker_sites$Name,y=marker_taxa$Taxa)

BLUPs_normalized = read_csv_retry(BLUP_normalized_path) %>%
  filter(Genotype %in% gt)
names(BLUPs_normalized)[c(1,4)]=c("X","BLUP")


# mxn => nxm
marker_geno= t(marker_geno)


if(trait == "RGB1_Plant_Avg_HEIGHT_MM"){
  dat = dat+1
}

BLUPS_foc_spec <- BLUPs_normalized %>%
  filter(Trait == trait) %>%
  filter(DAT == dat) %>%
  dplyr::select(X,BLUP)
BLUPS_foc_spec=BLUPS_foc_spec[order(BLUPS_foc_spec$X),]

CV_mat = read_csv_retry(CV_mat_file,rn = T)

full_df = data.frame()
nf=T




fold_vec=CV_mat[,run]
nfolds = max(fold_vec)
for(j in 1:nfolds){
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
                      DAT = dat,
                      run = run,
                      fold = j)
  
  if(nf){
    full_df = cur_df
    nf=F
  }else{
    full_df = rbind(full_df,cur_df)
  }
}


    





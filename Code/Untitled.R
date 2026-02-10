## ---------------------------
##
## Script name: GWAS_Shuffle_eval.R
##
## Purpose of script: Evaluating results of GWAS Shuffle
##
## Author: M.Sc. Robin Lindner
##
## Date Created: 2026-02-04
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

source("0_utils.R")   

## ---------------------------

## load up the packages we will need:  

input_dir = "../Data/Generated/GenomicPrediction/GWAS_Shuffle/"

files = list.files(input_dir)
full_df=data.frame()
nf=T
for(file in files){
  df=read.csv(paste0(input_dir,file),row.names = 1)
  if(nf){
    full_df = df
    nf=F
  }
  else{
    full_df = rbind(full_df,df)
  }
}


full_df$run = as.numeric(full_df$run)
full_df$fold = as.numeric(full_df$fold)
write.csv(full_df,paste0(input_dir,"full.csv"),row.names=F)


dat=15
trait = "RGB1_Plant_Avg_HEIGHT_MM"
subset <- full_df %>%
  filter(DAT==dat) %>%
  filter(trait == trait)
snp_tab = as.data.frame(table(subset$SNP_idxs))

res = GetIntersectingPartitionsForSNPSet(subset,c(snp_tab$Var1[1:2]))

GetIntersectingPartitionsForSNPSet <- function(subset,snps){
  validPartitions = data.frame()
  for(run in unique(subset$run)){
    for(fold in unique(subset$fold)){
      data_cut <- subset %>%
        filter(run == run & fold == fold)
      if(all(snps %in% data_cut$SNP_idxs)){
        cur = data.frame(Run = run,
                         Fold = fold)
        if(nf){
          validPartitions = cur
          nf=F
        }else{
          validPartitions = rbind(validPartitions,cur)
        }
      }
    }
  }
  return(validPartitions)
}
  
  
  
  

        
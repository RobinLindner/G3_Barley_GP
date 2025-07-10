## ---------------------------
##
## Script name: HPC_ESA
##
## Purpose of script: Extracting Significant Associations. Script to go through 
##                    results of 1456 GWAS and return the associations with 
##                    significance lower than set level.
##
## Author: M.Sc. Robin Lindner
##
## Date Created: 2025-07-08
##
## Copyright (c) Robin Lindner, 2025
## Email: robin.lindner@uni-potsdam.de
##
## ---------------------------
##
## Notes:
##   GWAS results need to follow the nomenclature: trait_._._dat
##
##   Genome-wide significance threshold was computed as 
##   0.05 / 1713 = 2.918856e-05
##
## ---------------------------

## set working directory for Mac and PC
# arg1 = folder of GWAS results
# arg2 = significance threshold
# arg3 = output file

args = commandArgs(trailingOnly = TRUE)
args = c("../Data/Generated/GWAS_results/",2.918856e-05,"../Data/Generated/") 
nf=T
sig_assoc=data.frame()
for(file in list.files(args[1])){
  trait = sub("_._._.+","",file)
  dat = sub(".FarmCPU.csv","",sub(".+_._._","",file))
  
  temp_result = read.csv(file.path(args[1],file))
  colnames(temp_result)[9] = "p_value"
  temp_result = temp_result %>%
    filter(p_value < as.numeric(args[2]))
  if(nrow(temp_result)>0){
    temp_result$DAT = dat
    temp_result$Trait = trait
    if(nf){
      sig_assoc = temp_result
      nf=F
    }else{
      sig_assoc = rbind(sig_assoc,temp_result)
    }
  }
}

write.csv(sig_assoc,args[3],row.names = F)


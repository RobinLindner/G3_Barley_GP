## ---------------------------
##
## Script name: HPC_GWAS
##
## Purpose of script: Execution of GWAS on a high performance cluster.
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
## Notes: This can be done locally or outsourced to handle high storage 
##        demand. 
##       
##       1456 measurement instances * 1.3Mb = 1.89 Gb
##   
##
## ---------------------------


## ---------------------------

## load up the packages we will need:  

library(rMVP)
library(bigmemory)
library(tidyr)
library(dplyr)
## ---------------------------

## load up our functions into memory:


## ---------------------------
# general time series HPC GWAS
# ---- cluster execution
# arg1 = input file of phenotypic data (normalized BLUPs)
# arg2 = output folder for GWAS results
# arg3 = vcf file of genotyping data
# arg4 = Kinship file
# arg5 = Supporting directory path

args = commandArgs(trailingOnly = TRUE)

## Read input data

BLUPs_normalized = read.csv(args[1])
phenotype = BLUPs_normalized %>%
  filter(!is.na(Value)) %>%
  pivot_wider(id_cols = Genotype,
              names_from = c(Trait,DAT),
              names_sep = "_._._",
              values_from = Value)

write.csv(phenotype,paste0(args[5],"/phenotype.csv"),row.names = F)

MVP.Data(fileVCF=args[3],filePhe = paste0(args[5],"/phenotype.csv"),fileKin = F, filePC = T,sep.phe = ",",out=paste0(args[5],"/B1K"))

K <-read.csv(args[4],row.names = 1)


datapref = paste0(args[5],"/B1K")

genotype  <- attach.big.matrix(paste0(datapref,".geno.desc"))     # Note: read the "<out>.imp" genotype
phenotype <- read.table(paste0(datapref,".phe"), header = TRUE)
map       <- read.table(paste0(datapref,".geno.map"), header = TRUE)
  

## Execute GWAS
for(i in 2:ncol(phenotype)){
  if(file.exists(paste0(args[2],"/",colnames(phenotype)[i],".FarmCPU.csv"))){
    print("skipped")
    next
  }
  print(colnames(phenotype)[i])
  MVP(
    phe=phenotype[, c(1, i)],
    geno=genotype,
    map=map,
    K=as.matrix(K),
    nPC.FarmCPU=4,
    vc.method="EMMA",
    maxLoop=10,
    method.bin="FaST-LMM",#"FaST-LMM","EMMA", "static"
    method=c("FarmCPU"),
    verbose = F,
    file.output=c("pmap"),
    outpath = args[2]
  )
}

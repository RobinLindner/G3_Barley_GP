## ---------------------------
##
## Script name: 0_utils.R
##
## Purpose of script: Contains all functions and file paths necessary for the 
##                    G3 Barley time-series GWAS publication.
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

## set working directory for Mac and PC

setwd(".")      # Robins's working directory (mac)


## ---------------------------

## load up the packages we will need:  
library(dplyr)
library(tidyr)
library(rMVP)
library(bigmemory)
library(car)
library(ggplot2)
library(ggpubr)
library(minpack.lm)
library(lme4)
library(bestNormalize)
library(lme4qtl)
library(caret)
library(Matrix)
library(MegaLMM)
library(tibble)

#library(topGO)
## ---------------------------
## load file paths into memory:

## !Need to be adjusted: 
# Depending on outpath of HPC_ESA.R call in 3_PhenomeWide_GWAS.R
sig_associations_file = "../Data/Generated/significant_associations.csv"

## Read-only paths:
vcf_file ="../Data/Genotype/B1K_final.vcf"
phenotype_nonHSR_file = "../Data/Phenotype/Merged_file_Tier1_Enviro.csv"
phenotype_HSR_file = "../Data/Phenotype/Spectrum_Tier1.csv"
trait_groups_file = "../Supplements/trait_groups.csv"

sites_file = "../Data/Genotype/B1K_final_sites.txt" # generated in TASSEL
taxa_file = "../Data/Genotype/B1K_final_taxa.txt" # generated in TASSEL
numeric_geno_file = "../Data/Genotype/B1K_final_numeric.txt" # generated in TASSEL

geno4GP_file = "../Supplements/GPgenotypes.txt"

ph_snp_ID_map = "../Supplements/ph_snp_map.csv"

## Read and write paths:
rMVP_out="../Data/Genotype/B1K_final"
geno_path = paste0(rMVP_out,".geno.desc")
map_path = paste0(rMVP_out,".geno.map")
geno_IDs_path = paste0(rMVP_out,".geno.ind")
GRM_path = "../Data/Genotype/B1K_final_GRM.csv"

BLUP_variance_components_path = "../Data/Generated/variance_components.csv"
BLUP_path = "../Data/Generated/BLUPs.csv"
BLUP_normalized_path = "../Data/Generated/BLUPs_normalized.csv"
phenotype_nonHSR_long_file = "../Supplements/Merged_file_Tier1_long.csv"

HSR_BLUP_variance_components_path = "../Data/Generated/HSR_variance_components.csv"
HSR_BLUP_path = "../Data/Generated/HSR_BLUPs.csv"
HSR_BLUP_normalized_path = "../Data/Generated/HSR_BLUPs_normalized.csv"
phenotype_HSR_long_file = "../Supplements/HSR_Merged_file_Tier1_long.csv"

geno_remap_file = "../Data/Genotype/B1K_SNP_remap.csv"

GP_CV_matrix_file = "../Supplements/GP_CV_mat2.csv"
GP_test_set_file = "../Supplements/GP_testSets.csv"
GP_valid_scenarios_file = "../Supplements/GP_valid_modelScenarios.csv"
## Write-only paths:
figure_dir = "../Figures/"
LD_out_file = "../Data/Genotype/LD/"
box_cox_parameter_file = "../Data/Generated/box_cox_parameters.csv"
HSR_box_cox_parameter_file = "../Data/Generated/HSR_box_cox_parameters.csv"
## ---------------------------
source("simpleM.R")
## load up functions into memory:

# For a selection of SNP IDs, return the candidate genes on the reference assembly
# @ param SNPs:       vector of SNP IDs
# @ param LD_ranges:  vector of chromosome-wise LD-decay length in bp
# @ param annotation: assembly annotation, mapping genes to positions on the genome
# @ param remap:      remapping of the SNP positions in the genotype date to a newer assembly.
GeneIDs_for_cand_SNP <- function(SNPs,LD_ranges,annotation,remap){
  snp_positions <- remap %>%
    filter(SNP %in% SNPs) %>%
    dplyr::select(new_position,Chromosome)
  if(length(SNPs)!= nrow(snp_positions)){
    print("At least one SNP not remapped!")
  }
  geneID_list = c()
  for(i in 1:nrow(snp_positions)){
    chrom = snp_positions$Chromosome[i]
    pos = snp_positions$new_position[i]
    qtl_start = pos - LD_ranges[chrom]
    qtl_end = pos + LD_ranges[chrom]
    chr_seqID = as.character(annotation$seqid[which(annotation$chromosome == paste0(chrom,"H"))])
    geneIDs = annotation %>%
      filter(start<= qtl_end & end >= qtl_start) %>%
      filter(seqid == chr_seqID) %>%
      dplyr::select(gene) %>%
      distinct(.keep_all = T) 
    geneID_list = c(geneID_list,geneIDs$gene)
    
  }
  return(geneID_list[-1])
}


GeneIDs_for_chromosome_range <- function(chrom,qtl_start,qtl_stop,annotation){
  geneID_list = c()
  chr_seqID = as.character(annotation$seqid[which(annotation$chromosome == paste0(chrom,"H"))])
  geneIDs = annotation %>%
    filter(start<= qtl_end & end >= qtl_start) %>%
    filter(seqid == chr_seqID) %>%
    dplyr::select(gene) %>%
    distinct(.keep_all = T) 
  geneID_list = c(geneID_list,geneIDs$gene)
  
  
  geneID_list=na.omit(unique(geneID_list))
  return(geneID_list)
}

FisherGOforAllDomains<-function(cand_genes,geneUniverse,geneID2GO_map, threshold){
  BP_GOdata <- new("topGOdata", ontology = "BP", allGenes = geneUniverse, annot = annFUN.gene2GO, gene2GO = geneID2GO_map)
  enriched_terms_BP = GetGO_pvalues(cand_genes = cand_genes,GOdata = BP_GOdata, domain="BP",threshold=threshold)
  
  MF_GOdata <- new("topGOdata", ontology = "MF", allGenes = geneUniverse, annot = annFUN.gene2GO, gene2GO = geneID2GO_map)
  enriched_terms_MF = GetGO_pvalues(cand_genes = cand_genes,GOdata = MF_GOdata, domain="MF",threshold=threshold)
  
  CC_GOdata <- new("topGOdata", ontology = "CC", allGenes = geneUniverse, annot = annFUN.gene2GO, gene2GO = geneID2GO_map)
  enriched_terms_CC = GetGO_pvalues(cand_genes = cand_genes,GOdata = CC_GOdata, domain="CC",threshold=threshold)
  
  return(rbind(enriched_terms_BP,enriched_terms_MF,enriched_terms_CC))
  
}

GetGO_pvalues <- function(cand_genes,GOdata,domain = "BP",threshold=0.05){
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  GOdata@allScores = factor(as.integer(geneNames %in% cand_genes))
  resultFisher_group <- getSigGroups(GOdata, test.stat)
  diagnostics_table = GenTable(object = GOdata,resultFisher_group,topNodes=min(1000,length(resultFisher_group@score)))
  diagnostics_table = diagnostics_table %>%
    filter(result1 < threshold)
  diagnostics_table$Domain = domain
  return(diagnostics_table)
}

prepareREVIGOIn <- function(GO_result,filename=""){
  tab = GO_result[,c(1,6)]
  write.table(tab,filename,sep = "\t",quote = F,row.names = F,col.names = F)
}  

GetTraitAtDAT_test <- function(single_tier_df,trait,dat){
  n_col=ncol(single_tier_df)
  trait_idx=which(colnames(single_tier_df)==trait)
  trait_df=single_tier_df[single_tier_df$DAT==dat,c(1:7,(n_col-15):n_col,trait_idx)]
  trait_df[,ncol(trait_df)]=as.numeric(trait_df[,ncol(trait_df)])
  trait_df=trait_df[!is.na(trait_df[,ncol(trait_df)]),]
  colnames(trait_df)[ncol(trait_df)] = "Trait"
  return(trait_df)
}

test_stdres_normality <- function(anova){
  # Extract the residuals
  aov_residuals <- residuals(object=anova)
  # Run Shapiro-Wilk test
  print(shapiro.test(x=aov_residuals))
}

pca_reduction <- function(matrix,threshold){
  na_rows = rowSums(is.na(matrix))>0
  print(paste0(sum(na_rows),"/",nrow(matrix)," rows contained NA values and were removed."))
  matrix = matrix[!na_rows,]
  pca =prcomp(matrix)
  exp_var = pca$sdev/sum(pca$sdev)
  cum_var = cumsum(exp_var)
  pc= min(which(cum_var>threshold))
  print(paste0(pc," PCs out of ",ncol(matrix)," variables, explain ",cum_var[pc]*100,"% of the variance"))
}
linearTemporalTrend <- function(data,tpCol,valueCol,groups=NA){
  if(!is.na(groups)){
    data = data %>% 
      filter(Group %in% groups)
  }
  x=as.matrix(cbind(rep(1,nrow(data)),data[,tpCol]))
  y=as.matrix(data[,valueCol])
  
  reg <- lm.fit(x,y)
  residuals = reg$residuals
  tss = sum((y-mean(y))^2)
  rss = sum(residuals^2)
  n = length(y)
  p = length(reg$coefficients)-1
  r_squared <- 1 - rss / tss
  adj_r_squared = 1 - (1 - r_squared) * ((n - 1) / (n - p - 1))
  if(!is.na(groups)){
    print(paste0("The linear trend for group(s): ",paste(groups,collapse=", ")," is ",reg$coefficients[2]," per day"))
  }else{
    print(paste0("The linear trend is ",reg$coefficients[2]," per day"))
    
  }
  print(paste0("Adjusted R^2 value: ", adj_r_squared))
}
linearTemporalTrend2 <- function(tp_vec,value_vec){
  
  x=as.matrix(cbind(rep(1,length(tp_vec)),tp_vec))
  y=as.matrix(value_vec)
  
  reg <- lm.fit(x,y)
  residuals = reg$residuals
  tss = sum((y-mean(y))^2)
  rss = sum(residuals^2)
  n = length(y)
  p = length(reg$coefficients)-1
  r_squared <- 1 - rss / tss
  adj_r_squared = 1 - (1 - r_squared) * ((n - 1) / (n - p - 1))
  
  print(paste0("The linear trend is ",reg$coefficients[2]," per day"))
  print(paste0("Adjusted R^2 value: ", adj_r_squared))
  return(list(reg$coefficients[1],reg$coefficients[2],adj_r_squared))
}

createTraitcomparisionDF <- function(snp,sig_assoc){
  geno <- read.table(numeric_geno_file)
  sites <- read.table(sites_file,header=T)
  taxa <- read.table(taxa_file,header=T)
  blups <- read.csv(BLUP_normalized_path)
  # find associated traits
  traits = sig_assoc %>%
    filter(SNP==snp) %>%
    dplyr::select(DAT,Trait) %>%
    distinct(.keep_all = T)
  
  # select the genotype of the SNP
  geno_at_snp=data.frame(Genotype=taxa$Taxa,alleleState=t(geno[which(sites$Name==snp),]))
  
  # filter out the unobserved genotypes
  filt_geno_at_snp = geno_at_snp %>%
    filter(Genotype %in% unique(blups$Genotype))
  colnames(filt_geno_at_snp)[2] ="AS"
  
  # filter the blups for measurement instances in which a significant association was detected
  relevant_blups <- blups %>%
    filter(Trait %in% traits$Trait & DAT %in% traits$DAT) %>%
    mutate(Trait=sub("FC1_Plant_","",Trait))
  
  # assign allele state as factor
  relevant_blups$Allele_state = as.factor(filt_geno_at_snp$AS[match(relevant_blups$Genotype,filt_geno_at_snp$Genotype)])
  
  stateNames = c("Homozygotic minor","Heterozygotic","Homozygotic major")
  
  relevant_blups <- relevant_blups %>%
    filter(!is.na(Allele_state)) %>%
    group_by(Trait,DAT) %>%
    mutate(NormValue = sym_min_max_scale(Value))
  
  
  relevant_blups$StateName="Homozygotic major"
  relevant_blups$StateName[relevant_blups$Allele_state==0] = stateNames[1]
  relevant_blups$StateName[relevant_blups$Allele_state==1] = stateNames[2]
  
  comparisons_list <- list(
    list(c(stateNames[0], stateNames[1])), 
    list(c(stateNames[1], stateNames[2])), 
    list(c(stateNames[0], stateNames[2]))
  )
  fake_x = seq(1,length(unique(relevant_blups$Trait)),1)
  
  
  ff=T
  for(trait in unique(relevant_blups$Trait)){
    test_blups <- relevant_blups %>%
      filter(Trait==trait)
    tukey = TukeyHSD(aov(NormValue~StateName,data=test_blups))
    temp = data.frame(Comparison=rownames(tukey$StateName),Diff =tukey$StateName[,1],p_value = tukey$StateName[,4])
    temp$Trait = trait
    rownames(temp)=NULL
    if(ff){
      stat_df = temp
      ff=F
    }else{
      stat_df = rbind(stat_df,temp)
    }
  }
  return(list(relevant_blups,stat_df))
}
sym_min_max_scale <- function(x){
  return(x/max(abs(x)))
}

jaccard <- function(set1,set2){
  return(length(intersect(set1,set2))/length(union(set1,set2)))
}

sameGroup <- Vectorize(function(Trait1,Trait2){
  trait_groups = read.csv(trait_groups_file)
  GroupsT1 = trait_groups$Group[match(Trait1,trait_groups$Trait)]
  GroupsT2 = trait_groups$Group[match(Trait2,trait_groups$Trait)]
  return(GroupsT1==GroupsT2)
})

#### ---- Genomic prediction ---- ####

UV_lm <- function(BLUPS_foc_spec,test_idx,H){
  
  K = as.matrix(read.csv(GRM_path,row.names = 1))
  
  train_data=BLUPS_foc_spec
  train_data$BLUP[test_idx]= NA
  
  test_y = BLUPS_foc_spec$BLUP[test_idx]
  
  test_geno = train_data$X[test_idx]
  train_geno = train_data$X[-test_idx]
  
  formula <- as.formula(BLUP~ (1|X))
  GBLUP_model <- relmatLmer(formula, train_data, relmat=list(X = K))
  
  pred_GBLUP_train <- predict(GBLUP_model)
  pred_GBLUP_test <- K[test_geno,train_geno] %*% solve(K[train_geno,train_geno]) %*% GBLUP_model@u
  pred_GBLUP_full=c()
  pred_GBLUP_full[na.omit(match(train_geno,BLUPS_foc_spec$X))]=pred_GBLUP_train
  pred_GBLUP_full[na.omit(match(test_geno,BLUPS_foc_spec$X))]=pred_GBLUP_test
  GBLUP_acc=cor(pred_GBLUP_test,test_y)
  
  print("GBLUP done")
  formula <- as.formula(BLUP~ (1|X))
  HBLUP_model <- relmatLmer(formula, train_data, relmat=list(X = H))
  
  pred_HBLUP_train <- predict(HBLUP_model)
  pred_HBLUP_test <- H[test_geno,train_geno] %*% solve(H[train_geno,train_geno]) %*% HBLUP_model@u
  pred_HBLUP_full=c()
  pred_HBLUP_full[na.omit(match(train_geno,BLUPS_foc_spec$X))]=pred_HBLUP_train
  pred_HBLUP_full[na.omit(match(test_geno,BLUPS_foc_spec$X))]=pred_HBLUP_test
  HBLUP_acc=cor(pred_HBLUP_test,test_y)
  
  print("HBLUP done")
  
  train_data_ext= cbind(train_data,data.frame(X_H = train_data$X))
  formula <- as.formula(BLUP~ (1|X)+(1|X_H))
  G_HBLUP_model <- relmatLmer(formula, train_data_ext, relmat=list(X = K,X_H=H))
  pred_GHBLUP_train <- predict(G_HBLUP_model)
  pred_GHBLUP_test <- K[test_geno,train_geno] %*% solve(K[train_geno,train_geno]) %*% G_HBLUP_model@u[1:length(train_geno)] + H[test_geno,train_geno] %*% solve(H[train_geno,train_geno]) %*% G_HBLUP_model@u[(length(train_geno)+1):length(G_HBLUP_model@u)]
  pred_GHBLUP_full=c()
  pred_GHBLUP_full[na.omit(match(train_geno,BLUPS_foc_spec$X))]=pred_GHBLUP_train[1:(length(pred_GHBLUP_train)/2)]
  pred_GHBLUP_full[na.omit(match(test_geno,BLUPS_foc_spec$X))]=pred_GHBLUP_test
  GHBLUP_acc=cor(pred_GHBLUP_test,test_y)
  
  print("GHBLUP done")
  
  acc=c(GBLUP_acc,HBLUP_acc,GHBLUP_acc)
  
  class = rep("Train",nrow(BLUPS_foc_spec))
  class[test_idx]="Test"
  cur_frame = data.frame(Geno = BLUPS_foc_spec$X,
                         BLUPs = BLUPS_foc_spec$BLUP,
                         GBLUP_Pred = pred_GBLUP_full,
                         HBLUP_Pred = pred_HBLUP_full,
                         GHBLUP_Pred = pred_GHBLUP_full,
                         Class = class)
  
  return(list(Accuracy=acc,Predictions=cur_frame))
}

UV_lm_MFE <- function(BLUPS_foc_spec,test_idx,snp_idxs,H,geno_mat){
  
  
  K = read.csv(GRM_path,row.names = 1)
  
  train_idx = setdiff(seq(1,215,1),test_idx)
  
  X_fe = NULL
  
  print(paste("Indices of significant SNP:",snp_idxs))
  geno = geno_mat[snp_idxs]
  
  gt = BLUPS_foc_spec$X
  
  # the Kinship matrix is created without the significant SNP to avoid redundancy
  Kinship_mat = MVP.K.VanRaden(as.big.matrix(as.matrix(geno_mat[gt,-snp_idxs])))
  rownames(Kinship_mat) = gt
  colnames(Kinship_mat) = gt
  
  # If the genotype of the significant SNPs has more than one position 
  # (should always be the case)
  if(ncol(geno)>0){
    X_fe = as.data.frame(geno[BLUPS_foc_spec$X[train_idx],])
    names(X_fe) = names(geno)
    if(ncol(geno)>1){
      red = reduceFE_matrix(X_fe) # If there is more than one SNP check for linear redundancy
      X_fe = red$mat
    }
    #fixed effect term is created as string
    fixed_effect_t = paste(paste(colnames(X_fe), collapse = " + "),"+")
    # the corresponding matrix is instantiated with NAs
    X_t = matrix(NA,nrow(BLUPS_foc_spec),ncol = ncol(X_fe))
    # The fixed effects for training genotypes are set 
    X_t[train_idx,]=as.matrix(X_fe)
    colnames(X_t) = names(X_fe)
    print(dim(X_t))
    print(dim(BLUPS_foc_spec))
    # training data is extended with the fixed effect columns for marker FE
    train_data=cbind(BLUPS_foc_spec,X_t)
  }
  else{
    fixed_effect_t =""
    train_data=BLUPS_foc_spec
  }
  
  # mask test genotypes in training dataset
  train_data$BLUP[test_idx]= NA
  
  # record test blup values for validation
  test_y = BLUPS_foc_spec$BLUP[test_idx]
  
  # get test and training genotypes
  test_geno = train_data$X[test_idx]
  train_geno = train_data$X[-test_idx]
  
  formula <- as.formula(paste0("BLUP ~ ", fixed_effect_t," (1|X)"))
  GBLUP_model <- relmatLmer(formula, train_data, relmat=list(X = Kinship_mat))
  
  pred_GBLUP_train <- predict(GBLUP_model)
  pred_GBLUP_test <- Kinship_mat[test_geno,train_geno] %*% solve(Kinship_mat[train_geno,train_geno]) %*% GBLUP_model@u
  pred_GBLUP_full=c()
  pred_GBLUP_full[na.omit(match(train_geno,BLUPS_foc_spec$X))]=pred_GBLUP_train
  pred_GBLUP_full[na.omit(match(test_geno,BLUPS_foc_spec$X))]=pred_GBLUP_test
  GBLUP_acc=cor(pred_GBLUP_test,test_y)
  
  print("GBLUP done")
  
  HBLUP_model <- relmatLmer(as.formula(paste0("BLUP ~ ", fixed_effect_t ," (1|X)")), train_data, relmat=list(X = H))
  
  pred_HBLUP_train <- predict(GBLUP_model)
  pred_HBLUP_test <- H[test_geno,train_geno] %*% solve(H[train_geno,train_geno]) %*% HBLUP_model@u
  pred_HBLUP_full=c()
  pred_HBLUP_full[na.omit(match(train_geno,BLUPS_foc_spec$X))]=pred_HBLUP_train
  pred_HBLUP_full[na.omit(match(test_geno,BLUPS_foc_spec$X))]=pred_HBLUP_test
  HBLUP_acc=cor(pred_HBLUP_test,test_y)
  
  print("HBLUP done")
  
  train_data_ext= cbind(train_data,data.frame(X_H = train_data$X))
  
  G_HBLUP_model <- relmatLmer(as.formula(paste0("BLUP ~ ", fixed_effect_t, " (1|X)+(1|X_H)")), train_data_ext, relmat=list(X = Kinship_mat,X_H=H))
  
  
  pred_GHBLUP_train <- predict(G_HBLUP_model)
  pred_GHBLUP_test <- Kinship_mat[test_geno,train_geno] %*% solve(Kinship_mat[train_geno,train_geno]) %*% G_HBLUP_model@u[1:length(train_geno)] + H[test_geno,train_geno] %*% solve(H[train_geno,train_geno]) %*% G_HBLUP_model@u[(length(train_geno)+1):length(G_HBLUP_model@u)]
  pred_GHBLUP_full=c()
  pred_GHBLUP_full[na.omit(match(train_geno,BLUPS_foc_spec$X))]=pred_GHBLUP_train[1:(length(pred_GHBLUP_train)/2)]
  pred_GHBLUP_full[na.omit(match(test_geno,BLUPS_foc_spec$X))]=pred_GHBLUP_test
  GHBLUP_acc=cor(pred_GHBLUP_test,test_y)
  
  print("GHBLUP done")
  GBLUP_acc_adj=GBLUP_acc
  HBLUP_acc_adj=HBLUP_acc
  GHBLUP_acc_adj=GHBLUP_acc
  if(!is.null(X_fe)){
    # if there are added fixed effects, adjust the prediction accuracy accordingly, 
    # by multiplying the marker state by its fixed effect value.
    test_X = as.data.frame(geno[BLUPS_foc_spec$X[test_idx],])
    
    if(ncol(test_X) != ncol(X_fe)){
      test_X = mapFullToRed(test_X,X_fe)
    }
    test_X = as.matrix(cbind(matrix(1,nrow(test_X),1),test_X))
    GBLUP_acc_adj=cor(pred_GBLUP_test + test_X %*% GBLUP_model@beta ,test_y)
    HBLUP_acc_adj=cor(pred_HBLUP_test + test_X %*% HBLUP_model@beta,test_y)
    GHBLUP_acc_adj=cor(pred_GHBLUP_test + test_X %*% G_HBLUP_model@beta,test_y)
  }
  
  if(!is.null(X_fe)){
    if(ncol(X_fe)>1){
      acc=c(GBLUP_acc_adj,HBLUP_acc_adj,GHBLUP_acc_adj,GBLUP_acc,HBLUP_acc,GHBLUP_acc,ncol(geno),red$n_red)
    }else{
      acc=c(GBLUP_acc_adj,HBLUP_acc_adj,GHBLUP_acc_adj,GBLUP_acc,HBLUP_acc,GHBLUP_acc,ncol(geno),0)
    }
  }else{
    acc=c(GBLUP_acc_adj,HBLUP_acc_adj,GHBLUP_acc_adj,GBLUP_acc,HBLUP_acc,GHBLUP_acc,ncol(geno),0)
  }
  
  class = rep("Train",nrow(BLUPS_foc_spec))
  class[test_idx]="Test"
  cur_frame = data.frame(Geno = BLUPS_foc_spec$X,
                         BLUPs = BLUPS_foc_spec$BLUP,
                         GBLUP_Pred = pred_GBLUP_full,
                         HBLUP_Pred = pred_HBLUP_full,
                         GHBLUP_Pred = pred_GHBLUP_full,
                         Class = class)
  
  if(!is.null(X_fe)){
    fe_frame_t = data.frame(FE=c("Intercept",colnames(X_fe)),
                            GBLUP_fe = GBLUP_model@beta,
                            HBLUP_fe = HBLUP_model@beta,
                            GHBLUP_fe = G_HBLUP_model@beta)
  }
  else{
    fe_frame_t = data.frame()
  }
  return(list(Accuracy=acc,Predictions=cur_frame,FE_sizes=fe_frame_t))
}






nFoldCV_lm_combined <- function(all_BLUPs,trait,dat,K,H,CV_mat,genotypes){
  Kinship_mat = as.matrix(K)
  H = as.matrix(H)
  gt = sort(genotypes)
  BLUPS_foc_spec <- all_BLUPs %>%
    filter(Trait == trait) %>%
    filter(DAT == dat) %>%
    filter(X %in% gt) %>%
    dplyr::select(X,BLUP)
  
  # reorder CV matrix & BLUPs (if necessary)
  CV_mat = CV_mat[match(gt,rownames(CV_mat)),]
  BLUPS_foc_spec=BLUPS_foc_spec[order(BLUPS_foc_spec$X),]
  
  result_frame = data.frame()
  acc_frame = data.frame(GBLUP_Accuracy=NA,HBLUP_Accuracy=NA,G_HBLUP_Accuracy=NA,Run=NA,Fold=NA)
  Kinship_mat = Kinship_mat[BLUPS_foc_spec$X,BLUPS_foc_spec$X]
  H = H[BLUPS_foc_spec$X,BLUPS_foc_spec$X]
  
  n_folds = length(unique(CV_mat[,1]))
  n_runs = ncol(CV_mat)
  for(i in 1:n_runs){
    fold_vec=CV_mat[,i]
    
    for(j in 1:n_folds){
      test_idx=which(fold_vec==j)
      
      train_data=BLUPS_foc_spec
      train_data$BLUP[test_idx]= NA
      
      test_y = BLUPS_foc_spec$BLUP[test_idx]
      
      test_geno = train_data$X[test_idx]
      train_geno = train_data$X[-test_idx]
      
      formula <- as.formula(BLUP~ (1|X))
      GBLUP_model <- relmatLmer(formula, train_data, relmat=list(X = Kinship_mat))
      
      pred_GBLUP_train <- predict(GBLUP_model)
      pred_GBLUP_test <- Kinship_mat[test_geno,train_geno] %*% solve(Kinship_mat[train_geno,train_geno]) %*% GBLUP_model@u
      pred_GBLUP_full=c()
      pred_GBLUP_full[na.omit(match(train_geno,BLUPS_foc_spec$X))]=pred_GBLUP_train
      pred_GBLUP_full[na.omit(match(test_geno,BLUPS_foc_spec$X))]=pred_GBLUP_test
      GBLUP_acc=cor(pred_GBLUP_test,test_y)
      
      print("GBLUP done")
      formula <- as.formula(BLUP~ (1|X))
      HBLUP_model <- relmatLmer(formula, train_data, relmat=list(X = H))
      
      pred_HBLUP_train <- predict(HBLUP_model)
      pred_HBLUP_test <- H[test_geno,train_geno] %*% solve(H[train_geno,train_geno]) %*% HBLUP_model@u
      pred_HBLUP_full=c()
      pred_HBLUP_full[na.omit(match(train_geno,BLUPS_foc_spec$X))]=pred_HBLUP_train
      pred_HBLUP_full[na.omit(match(test_geno,BLUPS_foc_spec$X))]=pred_HBLUP_test
      HBLUP_acc=cor(pred_HBLUP_test,test_y)
      
      print("HBLUP done")
      
      train_data_ext= cbind(train_data,data.frame(X_H = train_data$X))
      formula <- as.formula(BLUP~ (1|X)+(1|X_H))
      G_HBLUP_model <- relmatLmer(formula, train_data_ext, relmat=list(X = Kinship_mat,X_H=H))
      
      
      pred_GHBLUP_train <- predict(G_HBLUP_model)
      pred_GHBLUP_test <- Kinship_mat[test_geno,train_geno] %*% solve(Kinship_mat[train_geno,train_geno]) %*% G_HBLUP_model@u[1:length(train_geno)] + H[test_geno,train_geno] %*% solve(H[train_geno,train_geno]) %*% G_HBLUP_model@u[(length(train_geno)+1):length(G_HBLUP_model@u)]
      pred_GHBLUP_full=c()
      pred_GHBLUP_full[na.omit(match(train_geno,BLUPS_foc_spec$X))]=pred_GHBLUP_train[1:(length(pred_GHBLUP_train)/2)]
      pred_GHBLUP_full[na.omit(match(test_geno,BLUPS_foc_spec$X))]=pred_GHBLUP_test
      GHBLUP_acc=cor(pred_GHBLUP_test,test_y)
      
      print("GHBLUP done")
      
      acc_frame[(i-1)*n_folds+j,]=c(GBLUP_acc,HBLUP_acc,GHBLUP_acc,i,j)
      
      class = rep("Train",nrow(BLUPS_foc_spec))
      class[test_idx]="Test"
      cur_frame = data.frame(Geno = BLUPS_foc_spec$X,
                             BLUPs = BLUPS_foc_spec$BLUP,
                             GBLUP_Pred = pred_GBLUP_full,
                             HBLUP_Pred = pred_HBLUP_full,
                             GHBLUP_Pred = pred_GHBLUP_full,
                             Class = class,
                             Run = rep(i,nrow(BLUPS_foc_spec)),
                             Fold = rep(j,nrow(BLUPS_foc_spec)))
      if(i==1 && j==1){
        result_frame = cur_frame
      }else{
        result_frame = rbind(result_frame,cur_frame)
      }
    }
  }
  return(list(Accuracy=acc_frame,Predictions=result_frame))
}

partition <-function(idxs,folds){
  rem_idxs=idxs
  fold_vec=c()
  for(i in 1:(folds-1)){
    cur_idxs=sample(rem_idxs,length(idxs)/folds,replace=F)
    fold_vec[cur_idxs]=i
    rem_idxs=setdiff(rem_idxs,cur_idxs)
  }
  
  fold_vec[rem_idxs]=folds
  
  return(fold_vec)
}

reduceFE_matrix<- function(X){
  all_SNP = colnames(X)
  lc = findLinearCombos(X)
  all_FE = all_SNP
  for(comb in lc$linearCombos){
    all_FE[comb[2]] = paste0(all_FE[comb[2]],".",all_FE[comb[1]])
  }
  if(length(lc$remove)>0){
    all_FE = all_FE[-lc$remove]
    print(all_FE)
    X = X[,-lc$remove]
    print(X)
    colnames(X)=all_FE
  }
  return(list(mat=X,n_red = length(lc$remove)))
}

mapFullToRed <- function(full,red){
  i=1
  res_mat = matrix(NA,nrow=nrow(full),ncol=ncol(red))
  for(SNP_comb in names(red)){
    snps=strsplit(x = SNP_comb,split = ".",fixed = T)[[1]]
    res_mat[,i] = apply(as.data.frame(full[,snps]),1,mean)
    i=i+1
  }
  names(res_mat)=names(red)
  return(res_mat)
}

nFoldCV_lm_combined_MFE <- function(all_BLUPs,trait,dat,geno_mat,K,H,CV_mat,genotypes){
  
  Kinship_mat = as.matrix(K)
  H = as.matrix(H)
  gt = sort(genotypes)
  BLUPS_foc_spec <- all_BLUPs %>%
    filter(Trait == trait) %>%
    filter(DAT == dat) %>%
    filter(X %in% gt) %>%
    dplyr::select(X,BLUP)
  # reorder CV matrix & BLUPs (if necessary)
  CV_mat = CV_mat[match(gt,rownames(CV_mat)),]
  BLUPS_foc_spec=BLUPS_foc_spec[order(BLUPS_foc_spec$X),]
  
  result_frame = data.frame()
  fe_frame = data.frame()
  acc_frame = data.frame(GBLUP_Accuracy=NA,HBLUP_Accuracy=NA,G_HBLUP_Accuracy=NA,GBLUP_Accuracy_na=NA,HBLUP_Accuracy_na=NA,G_HBLUP_Accuracy_na=NA,Run=NA,Fold=NA,n_FE=NA,nr_redM = NA)
  Kinship_mat = Kinship_mat[BLUPS_foc_spec$X,BLUPS_foc_spec$X]
  H = H[BLUPS_foc_spec$X,BLUPS_foc_spec$X]
  
  n_folds = length(unique(CV_mat[,1]))
  n_runs = ncol(CV_mat)
  
  for(i in 1:n_runs){
    fold_vec=CV_mat[,i]
    
    for(j in 1:n_folds){
      test_idx=which(fold_vec==j)
      train_idx = which(fold_vec!=j)
      
      X_fe = NULL
      
      geno_train = geno_mat[BLUPS_foc_spec$X[train_idx],]
      K_train = MVP.K.VanRaden(as.big.matrix(as.matrix(geno_train)))
      #K_train = Kinship_mat[BLUPS_foc_spec$X[train_idx],BLUPS_foc_spec$X[train_idx]]
      rownames(K_train) = BLUPS_foc_spec$X[train_idx]
      colnames(K_train) = BLUPS_foc_spec$X[train_idx]
      
      Y_train = BLUPS_foc_spec[train_idx,]
      colnames(Y_train) = c("Taxa",trait)
      
      sig_snp = GetSignificantAssociationsForTrainingSet(Y_train = Y_train,
                                                          geno_train =  geno_train,
                                                          K_train = K_train)
      print(paste("Indices of significant SNP:",sig_snp))
      geno = geno_mat[sig_snp]
      
      
      Kinship_mat = MVP.K.VanRaden(as.big.matrix(as.matrix(geno_mat[gt,-sig_snp])))
      rownames(K_train) = gt
      colnames(K_train) = gt
      
      if(ncol(geno)>0){
        X_fe = as.data.frame(geno[BLUPS_foc_spec$X[train_idx],])
        names(X_fe) = names(geno)
        if(ncol(geno)>1){
          red = reduceFE_matrix(X_fe)
          X_fe = red$mat
        }
        fixed_effect_t = paste(paste(colnames(X_fe), collapse = " + "),"+")
        X_t = matrix(NA,nrow(BLUPS_foc_spec),ncol = ncol(X_fe))
        X_t[train_idx,]=as.matrix(X_fe)
        colnames(X_t) = names(X_fe)
        print(dim(X_t))
        print(dim(BLUPS_foc_spec))
        train_data=cbind(BLUPS_foc_spec,X_t)
      }
      else{
        fixed_effect_t =""
        train_data=BLUPS_foc_spec
      }
      
      train_data$BLUP[test_idx]= NA
      
      test_y = BLUPS_foc_spec$BLUP[test_idx]
      
      test_geno = train_data$X[test_idx]
      train_geno = train_data$X[-test_idx]
      
      formula <- as.formula(paste0("BLUP ~ ", fixed_effect_t," (1|X)"))
      GBLUP_model <- relmatLmer(formula, train_data, relmat=list(X = Kinship_mat))
      
      pred_GBLUP_train <- predict(GBLUP_model)
      pred_GBLUP_test <- Kinship_mat[test_geno,train_geno] %*% solve(Kinship_mat[train_geno,train_geno]) %*% GBLUP_model@u
      pred_GBLUP_full=c()
      pred_GBLUP_full[na.omit(match(train_geno,BLUPS_foc_spec$X))]=pred_GBLUP_train
      pred_GBLUP_full[na.omit(match(test_geno,BLUPS_foc_spec$X))]=pred_GBLUP_test
      GBLUP_acc=cor(pred_GBLUP_test,test_y)
      
      print("GBLUP done")
      
      HBLUP_model <- relmatLmer(as.formula(paste0("BLUP ~ ", fixed_effect_t ," (1|X)")), train_data, relmat=list(X = H))
      
      pred_HBLUP_train <- predict(GBLUP_model)
      pred_HBLUP_test <- H[test_geno,train_geno] %*% solve(H[train_geno,train_geno]) %*% HBLUP_model@u
      pred_HBLUP_full=c()
      pred_HBLUP_full[na.omit(match(train_geno,BLUPS_foc_spec$X))]=pred_HBLUP_train
      pred_HBLUP_full[na.omit(match(test_geno,BLUPS_foc_spec$X))]=pred_HBLUP_test
      HBLUP_acc=cor(pred_HBLUP_test,test_y)
      
      print("HBLUP done")
      
      train_data_ext= cbind(train_data,data.frame(X_H = train_data$X))
      
      G_HBLUP_model <- relmatLmer(as.formula(paste0("BLUP ~ ", fixed_effect_t, " (1|X)+(1|X_H)")), train_data_ext, relmat=list(X = Kinship_mat,X_H=H))
      
      
      pred_GHBLUP_train <- predict(G_HBLUP_model)
      pred_GHBLUP_test <- Kinship_mat[test_geno,train_geno] %*% solve(Kinship_mat[train_geno,train_geno]) %*% G_HBLUP_model@u[1:length(train_geno)] + H[test_geno,train_geno] %*% solve(H[train_geno,train_geno]) %*% G_HBLUP_model@u[(length(train_geno)+1):length(G_HBLUP_model@u)]
      pred_GHBLUP_full=c()
      pred_GHBLUP_full[na.omit(match(train_geno,BLUPS_foc_spec$X))]=pred_GHBLUP_train[1:(length(pred_GHBLUP_train)/2)]
      pred_GHBLUP_full[na.omit(match(test_geno,BLUPS_foc_spec$X))]=pred_GHBLUP_test
      GHBLUP_acc=cor(pred_GHBLUP_test,test_y)
      
      print("GHBLUP done")
      GBLUP_acc_adj=GBLUP_acc
      HBLUP_acc_adj=HBLUP_acc
      GHBLUP_acc_adj=GHBLUP_acc
      if(!is.null(X_fe)){
        test_X = as.data.frame(geno[BLUPS_foc_spec$X[test_idx],])
        
        if(ncol(test_X) != ncol(X_fe)){
          test_X = mapFullToRed(test_X,X_fe)
        }
        test_X = as.matrix(cbind(matrix(1,nrow(test_X),1),test_X))
        GBLUP_acc_adj=cor(pred_GBLUP_test + test_X %*% GBLUP_model@beta ,test_y)
        HBLUP_acc_adj=cor(pred_HBLUP_test + test_X %*% HBLUP_model@beta,test_y)
        GHBLUP_acc_adj=cor(pred_GHBLUP_test + test_X %*% G_HBLUP_model@beta,test_y)
      }
      
      if(!is.null(X_fe)){
        if(ncol(X_fe)>1){
          acc_frame[(i-1)*n_folds+j,]=c(GBLUP_acc_adj,HBLUP_acc_adj,GHBLUP_acc_adj,GBLUP_acc,HBLUP_acc,GHBLUP_acc,i,j,ncol(geno),red$n_red)
        }else{
          acc_frame[(i-1)*n_folds+j,]=c(GBLUP_acc_adj,HBLUP_acc_adj,GHBLUP_acc_adj,GBLUP_acc,HBLUP_acc,GHBLUP_acc,i,j,ncol(geno),0)
        }
      }else{
        acc_frame[(i-1)*n_folds+j,]=c(GBLUP_acc_adj,HBLUP_acc_adj,GHBLUP_acc_adj,GBLUP_acc,HBLUP_acc,GHBLUP_acc,i,j,ncol(geno),0)
      }
      
      class = rep("Train",nrow(BLUPS_foc_spec))
      class[test_idx]="Test"
      cur_frame = data.frame(Geno = BLUPS_foc_spec$X,
                             BLUPs = BLUPS_foc_spec$BLUP,
                             GBLUP_Pred = pred_GBLUP_full,
                             HBLUP_Pred = pred_HBLUP_full,
                             GHBLUP_Pred = pred_GHBLUP_full,
                             Class = class,
                             Run = rep(i,nrow(BLUPS_foc_spec)),
                             Fold = rep(j,nrow(BLUPS_foc_spec)))
      
      if(!is.null(X_fe)){
        fe_frame_t = data.frame(FE=c("Intercept",colnames(X_fe)),
                                GBLUP_fe = GBLUP_model@beta,
                                HBLUP_fe = HBLUP_model@beta,
                                GHBLUP_fe = G_HBLUP_model@beta,
                                Run = i,
                                Fold = j)
      }
      else{
        fe_frame_t = data.frame()
      }
      if(i==1 && j==1){
        result_frame = cur_frame
        fe_frame = fe_frame_t
      }else{
        result_frame = rbind(result_frame,cur_frame)
        fe_frame = rbind(fe_frame,fe_frame_t)
      }
    }
  }
  return(list(Accuracy=acc_frame,Predictions=result_frame,FE_sizes = fe_frame))
}


map_SNP_positions <- function(chr_map,chr_SNPs){
  result = data.frame(SNP = NA,old_position=NA,new_position=NA,Seq_Identity=NA)
  k=1
  for(i in 1:nrow(chr_SNPs)){
    old_pos = chr_SNPs$Position[i]
    selection <- chr_map %>%
      filter((S2<old_pos & E2>old_pos) | (S2>old_pos & E2<old_pos) ) 
    
    if(nrow(selection)==0){next}
    for(j in 1:nrow(selection)){
      per_base_comp_factor = selection$LEN1[j] / selection$LEN2[j]
      query_range=range(selection$S2[j],selection$E2[j])
      relative_pos = old_pos - query_range[1]
      new_pos = selection$S1[j] + relative_pos * per_base_comp_factor
      
      result[k,] = c(chr_SNPs$Name[i],old_pos,round(new_pos),selection$Identity...[j])
      k = k+1
      
    }
    
  }
  result$old_position = as.numeric(result$old_position)
  result$new_position = as.numeric(result$new_position)
  result$Seq_Identity = as.numeric(result$Seq_Identity)
  return(result)
}

#### ---- GWAS ---- ####
GetSignificantAssociationsForTrainingSet <- function(Y_train,geno_train,K_train){
 
  map = read.table(map_path,header = T)
  
  # compute number of PC that have more than 5% variance
  pca=prcomp(as.matrix(K_train))
  nPCs = tail(which((pca$sdev / sum(pca$sdev))>=0.05),1)
  
  # compute specific sig threshold
  genoPerChrom = SplitGenoPerChrom(as.matrix(geno_train),map)
  MeffTotal = SumMeffOverChrom(genoPerChrom)
  sprintf("Effective ratio: %f",MeffTotal/ncol(geno_train))
  sig_threshold = 0.05/MeffTotal
  
  # compute set-specific significant associations
  imMVP = MVP(
    phe=Y_train,
    geno=as.big.matrix(geno_train),
    map=map,
    K=as.matrix(K_train),
    nPC.FarmCPU=nPCs,
    vc.method="EMMA",
    maxLoop=10,
    method.bin="FaST-LMM",#"FaST-LMM","EMMA", "static"
    method=c("FarmCPU"),
    verbose = F,
    file.output = c()
  )
  p_values = imMVP$farmcpu.results[,3]
  print(paste("Number of significant SNP:",sum(p_values<=sig_threshold)))
  return(which(p_values<=sig_threshold))
}

SplitGenoPerChrom <- function(geno_mat,map){
  genoPerChrom = list()
  if(ncol(geno_mat)!=nrow(map)){
    print("Warning: Number of markers in the genotype matrix and in the map file are not equal!")
    print(ncol(geno_mat))
    print(nrow(map))
    return(genoPerChrom)
  }
  for(i in unique(map$CHROM)){
    genoPerChrom[[i]] = geno_mat[,map$CHROM==i]
  }
  return(genoPerChrom)
}

SumMeffOverChrom <- function(genoPerChrom){
  Meff=c()
  i=1
  for(genoMat in genoPerChrom){
    Meff[i]=computeMeff_mat(t(genoMat))
    i=i+1
  }
  return(sum(Meff))
}

CreateNewCV_mat <- function(filepath,nfold=10,nrun=100){
  gt = read.csv(geno4GP_file)$X
  CV_mat = replicate(nrun,partition(1:length(gt),nfold))
  rownames(CV_mat) = gt
  colnames(CV_mat) = paste0("run_",c(1:nrun))
  
  write.csv(CV_mat,filepath)
}

GetTestSetsForScenario <- function(scenarioID){
  temp = read.csv(GP_test_set_file,row.names = 1)
  split = t(data.frame(strsplit(temp[,paste0("X",scenarioID)],"_")))
  rownames(split) = NULL
  colnames(split) = c("Run","Fold")
  return(split)
}


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

read_table_retry <- function(file, max_retries = 3, delay = 5) {
  attempts <- 0
  while (attempts < max_retries) {
    attempts <- attempts + 1
    tryCatch({
      data <- read.table(file)
      return(data)
    }, error = function(e) {
      if (attempts == max_retries) stop("Max retries reached. Unable to read file.")
      Sys.sleep(delay) # Wait before retrying
    })
  }
}

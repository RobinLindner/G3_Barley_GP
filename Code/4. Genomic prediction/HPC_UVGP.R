setwd("Code")
source("0_utils.R")

## ---------------------------

## ---- create a directory for GP results ---- 

args = commandArgs(trailingOnly = T)

out_path = paste0("../",args[1])

trait = args[2]

dat = as.numeric(args[3])

CV_mat_file = paste0("../",args[4]) # contains all 100 splits

scenarioID = as.numeric(args[5]) # The ID of the scenario



# Read Kinship matrix as covariance structure K
K = read.csv(GRM_path,row.names = 1)

# Genotypes that had measurements for all time points & traits & HSR data
genotypes = read.table(geno4GP_file)$X

# Read numeric genotype to extract marker fixed effects (MFE)

marker_geno = t(as.matrix(read.table(numeric_geno_file_alt)))
marker_sites = read.table(sites_file_alt,header=T)
marker_taxa = read.table(taxa_file_alt,header=T)

rownames(marker_geno) = marker_taxa$Taxa
colnames(marker_geno) = marker_sites$Name
  
# Get the marker indices for the scenario
CV_mat = read.csv(GP_CV_matrix_file,row.names = 1)

validScenarios = read.csv(GP_valid_scenarios_file)
snp_idxs = as.numeric(unlist(strsplit(validScenarios$SNP_idxs[scenarioID],", ")))

# Read the 50 test sets from the supporting list
test_sets = GetTestSetsForScenario(scenarioID)

if(trait == "RGB1_Plant_Avg_HEIGHT_MM") {
  dat_foc = dat+1
}else{
  dat_foc = dat
}


# Prepare the normalized BLUEs (i.e. focal traits)
BLUEs_normalized = read.csv(BLUE_normalized_path) %>%
  filter(Genotype %in% genotypes) %>%
  filter(Trait == trait) %>%
  filter(DAT == dat_foc)

names(BLUEs_normalized)[c(1,4)]=c("X","BLUP")
BLUEs_normalized=BLUEs_normalized[match(rownames(CV_mat),BLUEs_normalized$X),]


# Prepare the normalized BLUPs of the HSR data (i.e. secondary traits)
HSR_BLUEs_normalized = read.csv(HSR_BLUE_normalized_path)%>%
  filter(Genotype %in% genotypes)



map <- read.table(map_path,header = T)

## ---- Univariate GP ----



for(i in 1:50){
  print(i)
  run = test_sets[i,1]
  fold= test_sets[i,2]
  
  X <- as.data.frame(marker_geno)
  colnames(X)=map$SNP
  
  dat_HS_blues <- HSR_BLUEs_normalized %>%
    filter(DAT==dat)
  colnames(dat_HS_blues)[1] = "X"
  print(colnames(dat_HS_blues))
  print(nrow(dat_HS_blues))
  HS_mat = pivot_wider(dat_HS_blues,id_cols = c(X),names_from = Trait,values_from = Value)%>%
    dplyr::select(where(~any(. !=1,na.rm=T))) %>%
    tibble::column_to_rownames(var="X")
  
  
  # prevent singularity
  rem=findLinearCombos(as.matrix(HS_mat))$remove
  HS_mat=HS_mat[-rem]
  
  HS_mat = HS_mat[rownames(CV_mat),]
  
  
  H=as.matrix(HS_mat) %*% t(as.matrix(HS_mat)) / ncol(HS_mat) 
  
  print(paste0("DAT:",dat," Trait:",trait))
  #print(paste0("number of sigificant SNPs: ",ncol(X)))


  fold_vec = CV_mat[,as.numeric(run)]
  test_idx = which(fold_vec==as.numeric(fold))

  
  # Returns a list containing 
  # $Accuracy:    | GBLUP Accuracy | HBLUP Accuracy | G+HBLUP Accuracy 
  # $Predictions: | Genotype | BLUPs | GBLUP Prediction | HBLUP Prediction | G+HBLUP Prediction | Test/Train
  res = UV_lm(BLUEs_normalized,test_idx,H)
  
  # Returns a list containing 
  # $Accuracy:    | GBLUP Accuracy | HBLUP Accuracy | G+HBLUP Accuracy | GBLUP Accuracy NonAdj | HBLUP Accuracy NonAdj | G+HBLUP Accuracy NonAdj | #Fixed effects | #reduced columns in genotype|
  # $Predictions: | Genotype | BLUPs | GBLUP Prediction | HBLUP Prediction | G+HBLUP Prediction | Test/Train
  # $FE_sizes:    | Fixed effect ID | GBLUP effect size | HBLUP effect size | G+HBLUP effect size 
  # NonAdj: prediction accuracy without considering marker effects cor(u,y) <=> cor(y',y)
  res_MFE = UV_lm_MFE(BLUEs_normalized,test_idx,snp_idxs,H,X)
  
  cur_acc = data.frame(GBLUP = res$Accuracy[1],
                       HBLUP = res$Accuracy[2],
                       GHBLUP = res$Accuracy[3],
                       Run = run,
                       Fold = fold,
                       Repetition = i)
  cur_pred = res$Predictions
  cur_pred$Run = run
  cur_pred$Fold = fold
  cur_pred$Repetition = i
  
  
  cur_acc_MFE = as.data.frame(t(data.frame(res_MFE$Accuracy)))
  colnames(cur_acc_MFE) = c("GBLUP","HBLUP","GHBLUP","GBLUPnonAdj","HBLUPnonAdj","GHBLUPnonAdj","nFixedEffects","nLinearlyDependentMarkers")
  cur_acc_MFE$Run = run
  cur_acc_MFE$Fold = fold
  cur_acc_MFE$Repetition = i
  
  cur_pred_MFE = res_MFE$Predictions
  cur_pred_MFE$Run = run
  cur_pred_MFE$Fold = fold
  cur_pred_MFE$Repetition = i
  
  cur_FE = res_MFE$FE_sizes
  cur_FE$Run = run
  cur_FE$Fold = fold
  cur_FE$Repetition = i
  
  if(i==1){
    scenario_acc_df = cur_acc
    scenario_pred_df = cur_pred
    scenario_MFE_acc_df = cur_acc_MFE
    scenario_MFE_pred_df = cur_pred_MFE
    scenario_MFE_fe_df = cur_FE

  }else{
    scenario_acc_df = rbind(scenario_acc_df,cur_acc)
    scenario_pred_df = rbind(scenario_pred_df,cur_pred)
    scenario_MFE_acc_df = rbind(scenario_MFE_acc_df,cur_acc_MFE)
    scenario_MFE_pred_df = rbind(scenario_MFE_pred_df,cur_pred_MFE)
    scenario_MFE_fe_df = rbind(scenario_MFE_fe_df,cur_FE)
    
  
  }
}
write.csv(scenario_acc_df,paste0(out_path,"/",trait,"_",dat,"_accuracy.csv"),row.names = F)
write.csv(scenario_pred_df,paste0(out_path,"/",trait,"_",dat,"_predictions.csv"),row.names = F)
write.csv(scenario_MFE_acc_df,paste0(out_path,"/",trait,"_",dat,"_MFE_accuracy.csv"),row.names = F)
write.csv(scenario_MFE_pred_df,paste0(out_path,"/",trait,"_",dat,"_MFE_predictions.csv"),row.names = F)
write.csv(scenario_MFE_fe_df,paste0(out_path,"/",trait,"_",dat,"_MFE_fixed_effects.csv"),row.names = F)
  
 
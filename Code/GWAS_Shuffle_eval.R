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
snp_tab = snp_tab[order(snp_tab$Freq,decreasing = T),]


rep_tab = as.data.frame(table(res$Run)) %>%
  filter(Freq==10)


res=GetIntersectingPartitionsForSNPSet(subset,1203)


t=Valid10foldCVrepetitions(subset,c(snp_tab$Var1[1]))

GetIntersectingPartitionsForSNPSet <- function(subset,snps){
  validPartitions = data.frame()
  for(Run in unique(subset$run)){
    for(Fold in unique(subset$fold)){
      data_cut <- subset %>%
        filter(run == Run) %>%
        filter(fold == Fold)
      if(all(snps %in% data_cut$SNP_idxs)){
        cur = data.frame(Run = Run,
                         Fold = Fold)
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

Valid10foldCVrepetitions <- function(subset,snps){
  res = GetIntersectingPartitionsForSNPSet(subset,snps)
  n_CVvalid = list()
  for(i in 1:10){
    valid_runs = as.data.frame(table(res$Run)) %>% 
      filter(Freq==i)
    n_CVvalid[[i]] = valid_runs$Var1
  }
  return(n_CVvalid)
}
traits = unique(full_df$Trait)
eval_df = data.frame()
nf=T
for(t in traits){
  trait_df = full_df %>%
    filter(Trait == t)
  for(dat in unique(trait_df$DAT)){
    cut_df = trait_df %>%
      filter(DAT == dat)
    snp_tab = as.data.frame(table(cut_df$SNP_idxs))
    snp_tab = snp_tab[order(snp_tab$Freq,decreasing = T),]
    for(i in 1:5){
      partitions = GetIntersectingPartitionsForSNPSet(cut_df,c(snp_tab$Var1[1:i]))
      CV10_part = Valid10foldCVrepetitions(cut_df,c(snp_tab$Var1[1:i]))
      print(length(CV10_part))
      if(length(CV10_part)!=0){
        
        splits = WhichSplitsFindSNPInAtLeastNFolds(CV10_part,5)
        cur = data.frame(Trait = t,
                         Dat = dat,
                         nSNP = i,
                         nPartitions = nrow(partitions),
                         nCV_valid = paste0(unlist(lapply(CV10_part,length)),collapse=", "),
                         SNP_idxs = paste0(snp_tab$Var1[1:i],collapse=", "),
                         CV10valid = paste0(CV10_part[[10]],collapse=", "),
                         AtLeast5folds = paste0(splits,collapse=", "))
      }else{
        cur = data.frame(Trait = t,
                         Dat = dat,
                         nSNP = i,
                         nPartitions = nrow(partitions),
                         nCV_valid = NA,
                         SNP_idxs = paste0(snp_tab$Var1[1:i],collapse=", "),
                         CV10valid = NA,
                         AtLeast5folds = NA)
      }
      
      if(nf){
        eval_df = cur
        nf=F
      }else{
        eval_df = rbind(eval_df,cur)
      }
    }
  }
}

WhichSplitsFindSNPInAtLeastNFolds <- function(Valid10FoldCVrep,n){
  reps = c()
  for(i in n:length(Valid10FoldCVrep)){
    reps = c(reps,Valid10FoldCVrep[[i]])
  }
  return(reps)
}
      

t = strsplit(eval_df$AtLeast5folds,", ")


run_df = data.frame()
for(i in 1:100){
  run = i
  cond = 0
  idxs = c()
  for(j in 1:length(t)){
    if(run %in% t[[j]]){
      cond = cond + 1
      idxs = c(idxs,j)
    }
  }
  
  scenarios = eval_df[idxs,1:3]
  scenarios$run = i
  scenarios$ncond = cond
  if(nf){
    run_df = scenarios
    nf=F
  }else{
    run_df = rbind(run_df,scenarios)
  }
}

for(i in 1:100){
  run_subdf = run_df %>%
    filter(run == i)
  print(paste0("Run ",i,": ",length(unique(paste0(run_subdf$Trait,"_",run_subdf$Dat)))))
  if(length(unique(paste0(run_subdf$Trait,"_",run_subdf$Dat)))==15){
    sprintf("Split %d finds at least one SNP in at least 5 training sets in all conditions.",i)
  }
}


possible_split = run_df %>%
  filter(run==14)


# Given a number of valid training sets that find the same associations (n/10)
#   1. find which splits find those associations in at least the given 
#     number of training sets.
#   2. Of those, rank which split fullfill this condition for the most scenarios (Trait_DAT) 
#   3. return the scenarios of that split, with the indices of the test sets omissions
#     that find the associations.
#   
#   This will show which splits are the most promising for any number of valid training sets we choose.
#   Additionally this will show if the training sets that find the SNP associations are shared between scenarios,
#   Potentially allowing us to choose the same training sets across scenarios-

FindBestRun<- function(full_df,nFolds){
  
  ## eval_df: | Trait | DAT | nSNP | nPartitions | nCV_valid | SNP_idx | CV10valid | AtLeastNFolds
  ## 
  ## For each scenario (trait at dat) we include up to 5 most-frequent SNP
  ## and record how many of the 1000 training sets found each set of SNP as 
  ## significantly associated.
  ## 
  ## nCV_valid: is a 10 element vector recording how many of the 100 splits are 
  ##            valid for a n-foldCV *
  ## CV10valid: is a comma seperated list containing the splits that are valid 
  ##            for a 10 fold CV Ü
  ## AtLeastNFolds: is a comma seperated list containing the splits that are 
  ##            valid for at least a n-fold CV (user input) *
  ## 
  ## * For the a specific SNP-set and scenario

  
  traits = unique(full_df$Trait)
  eval_df = data.frame()
  nf=T
  for(t in traits){
    trait_df = full_df %>%
      filter(Trait == t)
    for(dat in unique(trait_df$DAT)){
      cut_df = trait_df %>%
        filter(DAT == dat)
      snp_tab = as.data.frame(table(cut_df$SNP_idxs))
      snp_tab = snp_tab[order(snp_tab$Freq,decreasing = T),]
      for(i in 1:5){
        partitions = GetIntersectingPartitionsForSNPSet(cut_df,c(snp_tab$Var1[1:i]))
        CV10_part = Valid10foldCVrepetitions(cut_df,c(snp_tab$Var1[1:i]))
        if(length(CV10_part)!=0){
          
          splits = WhichSplitsFindSNPInAtLeastNFolds(CV10_part,nFolds)
          cur = data.frame(Trait = t,
                           Dat = dat,
                           nSNP = i,
                           nPartitions = nrow(partitions),
                           nCV_valid = paste0(unlist(lapply(CV10_part,length)),collapse=", "),
                           SNP_idxs = paste0(snp_tab$Var1[1:i],collapse=", "),
                           CV10valid = paste0(CV10_part[[10]],collapse=", "),
                           AtLeastNfolds = paste0(splits,collapse=", "))
        }else{
          cur = data.frame(Trait = t,
                           Dat = dat,
                           nSNP = i,
                           nPartitions = nrow(partitions),
                           nCV_valid = NA,
                           SNP_idxs = paste0(snp_tab$Var1[1:i],collapse=", "),
                           CV10valid = NA,
                           AtLeastNfolds = NA)
        }
        
        if(nf){
          eval_df = cur
          nf=F
        }else{
          eval_df = rbind(eval_df,cur)
        }
      }
    }
  }
  
  # Given the first evaluation we iterate through the runs and get the number of scenarios 
  
  #  i.e. the same split leads to significant SNP at all traits & time points
  # conditional on the split leading to at least N training sets that find the 
  # SNP-set.
  
  # returns Run_df:
  # |Run| nTnD | maxSNPperTnD
  #
  # nTnD: number of trait-timepoint combinations where the split lead to at 
  #       least n valid training sets 
  # v: 15 element vector containing the max number of SNP for the TnD
  
  splitsWithNFoldsInAllScenarios = strsplit(eval_df$AtLeastNfolds,", ")
  
  traits=c(rep("RGB1_Plant_Avg_HEIGHT_MM",5),
          rep("SC_Plant_Weight",5),
          rep("VNIR_Plant_NDVI.avg",5))
  dats = c(seq(15,43,7),
          seq(14,42,7),
          seq(14,42,7))
  scenario_lookup = paste0(traits,"_",dats)
  
  curr_best = 10
  best_split = data.frame()
  run_df = data.frame()
  for(i in 1:100){
    run = i
    idxs = c()
    for(j in 1:length(splitsWithNFoldsInAllScenarios)){
      if(run %in% splitsWithNFoldsInAllScenarios[[j]]){
        idxs = c(idxs,j)
      }
    }
    
    scenarios = eval_df[idxs,1:3]
    grouped_by_TnD=scenarios %>% dplyr::group_by(Trait,Dat) %>% summarize(maxSNPperDat = max(nSNP))
    TnD = paste0(grouped_by_TnD$Trait,"_",grouped_by_TnD$Dat)
    SNPperTnD = c()
    SNPperTnD[match(TnD,scenario_lookup)] = grouped_by_TnD$maxSNPperDat
    
    if(nrow(grouped_by_TnD)>curr_best){
      best_split = eval_df[idxs,]
      curr_best = nrow(grouped_by_TnD)
      best_split_ID = i
    }
    
    if(nrow(grouped_by_TnD)){
      sprintf("Split %d finds at least one SNP in at least %d training sets in all conditions.",run,nFolds)
    }
    
    cur = data.frame(Run = run,
                     nTnD = nrow(grouped_by_TnD),
                     maxSNPperTnD = paste(SNPperTnD,collapse=", ")
                     )
    if(nf){
      run_df = cur
      nf=F
    }else{
      run_df = rbind(run_df,cur)
    }
  }
  
  
  # Give more information about the data splits that are valid for a maximal
  # number scenarios
  # best_split = #ID
  testing_folds = c()
  for(i in 1:nrow(best_split)){
    subset = full_df %>%
      filter(run == best_split_ID) %>%
      filter(Trait == best_split$Trait[i]) %>%
      filter(DAT == best_split$Dat[i]) %>%
      filter(SNP_idxs %in% unlist(strsplit(best_split$SNP_idxs[i],", ")))
    folds = paste0(unique(subset$fold), collapse = ", ")
    testing_folds[i] = folds
  }
  best_split$testingFolds=testing_folds
  drop_columns = c("nPartitions","nCV_valid","CV10valid","AtLeastNfolds")
  best_split = best_split[,!colnames(best_split) %in% drop_columns]
  
  
  return(list(best_split_ID,eval_df,run_df,best_split))
  
}

WithMin5fold = FindBestRun(full_df,5)
WithMin4fold = FindBestRun(full_df,4)
WithMin3fold = FindBestRun(full_df,3)

WithMin1fold = FindBestRun(full_df,1)
WithMin2fold = FindBestRun(full_df,2)

eval_df_3fold= WithMin3fold[[2]]
run_df_3fold = WithMin3fold[[3]]
best_run_3fold = WithMin3fold[[4]]


GetSpecRunDf <- function(run_ID,full_df,eval_df){
  splitsWithNFoldsInAllScenarios = strsplit(eval_df$AtLeastNfolds,", ")
  idxs=c()
  for(j in 1:length(splitsWithNFoldsInAllScenarios)){
    if(run_ID %in% splitsWithNFoldsInAllScenarios[[j]]){
      idxs = c(idxs,j)
    }
  }
  srun_df = eval_df[idxs,]
  testing_folds = c()
  for(i in 1:nrow(srun_df)){
    subset = full_df %>%
      filter(run == run_ID) %>%
      filter(Trait == srun_df$Trait[i]) %>%
      filter(DAT == srun_df$Dat[i]) %>%
      group_by(Trait,DAT,fold) %>%
      summarise(check = all(unlist(strsplit(srun_df$SNP_idxs[i],", ")) %in% SNP_idxs))
    
    folds = paste0(unique(subset$fold[subset$check]), collapse = ", ")
    testing_folds[i] = folds
  }
  srun_df$testingFolds=testing_folds
  drop_columns = c("nPartitions","nCV_valid","CV10valid","AtLeastNfolds")
  srun_df = srun_df[,!colnames(srun_df) %in% drop_columns]
}


run1_3fold = GetSpecRunDf(2,full_df,eval_df_3fold)

# 20 x 3/10 fold-CV:
# 20 splits were selected that each contained 3 training sets that after 
# independent GWAS, lead to inclusion of the same marker fixed effects for all
# traits and time points. These training sets are not the same between scenarios.


# For 4fold-CV only 2 data splits fullfill this condition ->
# Using a higher number of folds (>3), one scenario, namely pot weight at 35 DAT, 
# would have to be excluded or scenarios would have to be evaluated with
# different numbers of training-testing splits.

run_df_1fold = WithMin1fold[[3]]
eval_df_1fold = WithMin1fold[[2]]

run_df_2fold = WithMin2fold[[3]]
eval_df_2fold = WithMin2fold[[2]]

run1_1fold = GetSpecRunDf(1,full_df,eval_df_1fold)


best_test_sets = GetBestTestSets(full_df,eval_df_1fold,run_df_1fold)

GetBestTestSets <- function(full_df,eval_df,run_df){
  final_set = list()
  traits=c(rep("RGB1_Plant_Avg_HEIGHT_MM",5),
           rep("SC_Plant_Weight",5),
           rep("VNIR_Plant_NDVI.avg",5))
  dats = c(seq(15,43,7),
           seq(14,42,7),
           seq(14,42,7))
  scenario_lookup = paste0(traits,"_",dats)
  
  for(scenario in scenario_lookup){
    
  }
  for(i in 1:nrow(run_df)){
    # select only ones valid for all scenarios
    if(as.numeric(run_df$nTnD[i])<15){
      next
    }
    spec_run_df = GetSpecRunDf(i,full_df,eval_df)
    t = strsplit(spec_run_df$testingFolds,", ")
    valid_partitions = c()
    for(j in 1:10){
      subset = spec_run_df[unlist(lapply(t,function(x) j %in% x)),]
      if(!all(scenario_lookup %in% paste0(subset$Trait,"_",subset$DAT))){
        next
      }else{
        valid_partitions = c(valid_partitions,j)
      }
      
    }
    final_set[[i]] = valid_partitions
    
  }
}


usedSNP = eval_df_1fold %>%
  filter(!is.na(nCV_valid)) %>%
  dplyr::select(Trait,Dat,nSNP,SNP_idxs,AtLeastNfolds) %>%
  distinct(.keep_all = T)

usedSNP$nAtLeastNfolds = unlist(lapply(strsplit(usedSNP$AtLeastNfolds, ", "),FUN=length))

validScenarios = usedSNP %>%
  filter(nAtLeastNfolds>=50) %>%
  group_by(Trait,Dat) %>%
  filter(nSNP == max(nSNP))

write.csv(validScenarios,GP_valid_scenarios_file,row.names = F)


validTestSets = data.frame()
nf=T
for(i in 1:100){ # iterate through all runs
  for(j in 1:10){ # iterate through all folds
    subset <- full_df %>%
      filter(run == i) %>%
      filter(fold == j)
    
    scenarios =c()
    for(k in 1:nrow(validScenarios)){ # if the SNPs for the scenario are predicted, add the scenario 
      subsubset <- subset %>%
        filter(Trait == validScenarios$Trait[k]) %>%
        filter(DAT == validScenarios$Dat[k])
      snp_idxs = unlist(strsplit(validScenarios$SNP_idxs[k],", "))
      if(all(snp_idxs %in% subsubset$SNP_idxs)){
        scenarios[k]=1
      }else{
        scenarios[k]=0
      }
    
    }
    if(sum(scenarios)==0){next}
    cur=data.frame(run = i,fold=j,
                   nScenarios = sum(scenarios),
                   scenarios = paste0(scenarios,collapse = ", "))
    if(nf){
      validTestSets = cur
      nf=F
    }else{
      validTestSets = rbind(validTestSets,
                           cur)
    }
  }
}

validTestSets=validTestSets[order(validTestSets$nScenarios,decreasing = T),]

testSetsForScenario <- vector(mode = "list", length = 15)
TS_count_perScenario = rep(0,15)
i=0
while(any(TS_count_perScenario<50)){
  i=i+1
  t = unlist(strsplit(validTestSets$scenarios[i],", "))
  for(j in 1:length(t)){
    if(t[j] == 1 & TS_count_perScenario[j]<50){
      testSetsForScenario[[j]] = c(testSetsForScenario[[j]], paste0(validTestSets$run[i],"_",validTestSets$fold[i]))
      TS_count_perScenario[j] = TS_count_perScenario[j] + 1
    }
  }
}

GetTestSetsForScenario <- function(scenarioID){
  temp = read.csv(GP_test_set_file,row.names = 1)
  split = t(data.frame(strsplit(temp[,paste0("X",scenarioID)],"_")))
  rownames(split) = NULL
  colnames(split) = c("Run","Fold")
  return(split)
}


t=data.frame(testSetsForScenario)
colnames(t) = seq(1,15,1)
write.csv(t,file = "../Supplements/GP_testSets.csv")

GetTestSetsForScenario(1)


#### plot genotypes used in training
CV_mat = read.csv(GP_CV_matrix_file,row.names = 1)
sets = GetTestSetsForScenario(1)
row = 1
run = as.numeric(sets[row,1])
fold = as.numeric(sets[row,2])
genotypes =rownames(CV_mat)[CV_mat[,run] != fold]

get_genotypes_used_in_training = function(scenarioID){
  sets = GetTestSetsForScenario(scenarioID)
  genotypes = c()
  for(row in 1:nrow(sets)){
    run = as.numeric(sets[row,1])
    fold = as.numeric(sets[row,2])
    genotypes =c(genotypes,rownames(CV_mat)[CV_mat[,run] != fold])
  }
  frequency_df=as.data.frame(table(genotypes))
  colnames(frequency_df) = c("Genotype",paste0("sc_",scenarioID))
  return(frequency_df)
}

for(i in 1:15){
  df = get_genotypes_used_in_training(i)
  if(i==1){
    full_df = df
  }
  else{
    full_df = merge(full_df,df)
  }
}
full_df_wo_PW = full_df[,c(2:6,11:16)]

rm = rowMeans(full_df_wo_PW[,-1])

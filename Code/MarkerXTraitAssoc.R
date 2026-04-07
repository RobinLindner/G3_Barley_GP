## ---------------------------
##
## Script name: 
##
## Purpose of script: Marker X Trait interactions i.e. pleiotropy testing.
##
## Author: M.Sc. Robin Lindner
##
## Date Created: 2026-02-18
##
## Copyright (c) Robin Lindner, 2026
## Email: robin.lindner@uni-potsdam.de
##
## ---------------------------
##
## Notes:
#Significant Marker effect, Non-sig Interaction	
# Pleiotropy: 
# The marker affects both traits in a similar way (common genetic control).

# Significant Marker X Trait Interaction
# Differential Effect: 
# The marker has a significantly stronger (or opposite) effect on one trait compared to the other.
##   
##
## ---------------------------

## set working directory for Mac

source("0_utils.R") 

## For pSNP
pSNP1 = "chr4H_632274504"
pSNP2 = "chr2H_704876309"
pSNP3 = "chr1H_349277204"


# load data
MVP.Data(fileVCF = "../Data/Genotype/B1K_red.vcf",out = "B1K_red")

geno = as.matrix(attach.big.matrix("../Data/Genotype/B1K_red.geno.desc"))
map = read.table(file = "../Data/Genotype/B1K_red.geno.map",header=T)
ind = read.table(file = "../Data/Genotype/B1K_red.geno.ind",header=F)
rownames(geno) = transformTaxaToGenotypeID(ind$V1)
colnames(geno) = map$SNP

K_red = MVP.K.VanRaden(as.big.matrix(geno))
colnames(K_red) = rownames(geno)
rownames(K_red) = rownames(geno)
write.csv(K_red,"../Data/Genotype/B1K_GRM_red.csv")

K_red = read.csv("../Data/Genotype/B1K_GRM_red.csv",row.names = 1)

normalized_BLUEs = read.csv(BLUE_normalized_path)

assoc_frame = read.csv("../Data/Generated/significant_associations.csv")

# 1. Get associated traits for pSNP
pSNP = pSNP2


SolveMarkerXTrait = function(pSNP){
  pSNP_trait = assoc_frame %>%
    filter(SNP == pSNP) %>%
    dplyr::select(Trait,DAT) %>%
    distinct(.keep_all = T) %>%
    mutate(TnD = paste(Trait,DAT,sep="_"))
  
  
  # 2. Stack BLUPs of these traits 
  relevant_BLUE = normalized_BLUEs %>%
        filter(Trait %in% pSNP_trait$Trait)
  
  # 3. Append Marker state 
  marker_vec = geno[,pSNP]
  relevant_BLUE$MarkerState = marker_vec[match(relevant_BLUE$Genotype,rownames(geno))]
  relevant_BLUE=relevant_BLUE%>%
    filter(!is.na(MarkerState)) %>%
    filter(!Trait %in% c("SC_Plant_Weight"))
  
  
  #4. Build model for MarkerXTrait
  # 4.1 separate models for different DAT
  dat = unique(relevant_BLUE$DAT)
  allCoeff = data.frame()
  nf=T
  for(d in dat){
    tp_rel_BLUE = relevant_BLUE%>%filter(DAT==d)
    if(!length(unique(tp_rel_BLUE$Trait))>1){
      print(paste0("TP: ",d," not eligible for model"))
    }
    else{
      model = lm(Value~MarkerState*Trait,data = tp_rel_BLUE)
      cur = as.data.frame(anova(model))
      cur$eta = cur$`Sum Sq` / (cur$`Sum Sq` + cur$`Sum Sq`[4])
      cur$DAT = d
      cur$Effect = rownames(cur)
      if(nf){
        allCoeff = cur
        nf=F
      }else{
        allCoeff = rbind(allCoeff,
                         cur)
      }
    }
  }
  rownames(allCoeff) = NULL
  allCoeff=allCoeff[order(allCoeff$DAT,decreasing = F),]
  plot_df = allCoeff[allCoeff$Effect!="Residuals",]
  p<-ggplot(plot_df,aes(x=DAT,y=eta,color=Effect))+
    geom_line()
  print(p)
  nf=T
  for(dat in sort(unique(allCoeff$DAT))){
    traits = unique(relevant_BLUE$Trait[relevant_BLUE$DAT == dat])
    cur = data.frame(DAT = dat,
                     eta_MxT = allCoeff$eta[allCoeff$DAT==dat & allCoeff$Effect=="MarkerState:Trait"],
                     eta_M = allCoeff$eta[allCoeff$DAT==dat & allCoeff$Effect=="MarkerState"],
                     eta_T = allCoeff$eta[allCoeff$DAT==dat & allCoeff$Effect=="Trait"],
                     p_MxT = allCoeff$`Pr(>F)`[allCoeff$DAT==dat & allCoeff$Effect=="MarkerState:Trait"],
                     p_M = allCoeff$`Pr(>F)`[allCoeff$DAT==dat & allCoeff$Effect=="MarkerState"],
                     p_T = allCoeff$`Pr(>F)`[allCoeff$DAT==dat & allCoeff$Effect=="Trait"],
                     nTrait = length(traits),
                     Trait = paste(traits,collapse = ", "))
    if(nf){
      full_df = cur
      nf=F
    }
    else{
      full_df = rbind(full_df,
                      cur)
    }
  }
  print(sprintf("%d models were evaluated.",length(unique(full_df$DAT))))
  print(sprintf("%d models are significantly affected by the allele state",sum(full_df$p_M<0.05)))
  print(sprintf("%d models are significantly affected by the interaction term",sum(full_df$p_MxT<0.05)))
  print(sprintf("Allele state eta squared higher in %d models",sum(full_df$eta_M>full_df$eta_MxT)))
  print(sprintf("Interaction eta squared higher in %d models",sum(full_df$eta_MxT>full_df$eta_M)))
}

SolveMarkerXTrait(pSNP1)
SolveMarkerXTrait(pSNP2)
SolveMarkerXTrait(pSNP3)


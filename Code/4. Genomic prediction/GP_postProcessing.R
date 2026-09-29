## ---------------------------
##
## Script name: GP_eval_2
##
## Purpose of script: Summarize and visualize the new GP results
##
## Author: M.Sc. Robin Lindner
##
## Date Created: 2026-02-11
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

# load Uni-variate results 
files = list.files("../Data/Generated/GenomicPrediction/UVGP_old//")

UV_accuracy_df = data.frame()
UV_prediction_df = data.frame()
UV_FE_df = data.frame()

acc_nf = T
pred_nf= T
fe_nf = T
for(file in files){
  if(grepl("_accuracy.csv",file)){
    if(grepl("_MFE_",file)){
      scenario = sub("_MFE_accuracy.csv","",file)
    }else{
      scenario = sub("_accuracy.csv","",file)
    }
    t = strsplit(scenario,"_")
    dat = t[[1]][length(t[[1]])]
    trait = paste0((t[[1]][-length(t[[1]])]),collapse = "_")
    df = read.csv(paste0("../Data/Generated/GenomicPrediction/UVGP_old//",file))
    df$Trait = trait
    df$Dat = dat
    df$MFE = grepl("_MFE_",file)
    if(acc_nf){
      UV_accuracy_df = df
      acc_nf = F
    }else{
      UV_accuracy_df = bind_rows(UV_accuracy_df, df)
    }
  }
  else if(grepl("_predictions.csv", file)){
    if(grepl("_MFE_", file)){
      scenario = sub("_MFE_predictions.csv","",file)
    }else{
      scenario = sub("_predictions.csv","",file)
    }
    t = strsplit(scenario,"_")
    dat = t[[1]][length(t[[1]])]
    trait = paste0((t[[1]][-length(t[[1]])]),collapse = "_")
    df = read.csv(paste0("../Data/Generated/GenomicPrediction/UVGP/",file))
    df$Trait = trait
    df$Dat = dat
    df$MFE = grepl("_MFE_", file)
    if(pred_nf){
      UV_prediction_df = df
      pred_nf = F
    }else{
      UV_prediction_df = bind_rows(UV_prediction_df, df)
    }
  }
  else if(grepl("_fixed_effects", file)){
    scenario = sub("_MFE_fixed_effects.csv","",file)
    t = strsplit(scenario,"_")
    dat = t[[1]][length(t[[1]])]
    trait = paste0((t[[1]][-length(t[[1]])]),collapse = "_")
    df = read.csv(paste0("../Data/Generated/GenomicPrediction/UVGP/",file))
    df$Trait = trait
    df$Dat = dat
    if(fe_nf){
      UV_FE_df = df
      fe_nf = F
    }else{
      UV_FE_df = bind_rows(UV_FE_df, df)
    }
  }
}

prelimAcc <- UV_accuracy_df %>%
  dplyr::select(GBLUP,HBLUP,GHBLUP,Trait,Dat,MFE) %>%
  group_by(Trait,Dat,MFE) %>%
  summarize(meanGBLUP = mean(GBLUP,na.rm=T),
            meanHBLUP = mean(HBLUP,na.rm=T),
            meanGHBLUP = mean(GHBLUP,na.rm=T))

plotDF <- UV_accuracy_df %>%
  dplyr::select(GBLUP,HBLUP,GHBLUP,Trait,Dat,MFE) %>%
  pivot_longer(cols = c(GBLUP,HBLUP,GHBLUP),names_to = "Model",values_to = "PredictionAccuracy")


h2_df = read.csv("../Data/Generated/nonHSR_h2.csv",row.names = 1)

scenarios = c(paste(seq(15,43,7),c("RGB1_Plant_Avg_HEIGHT_MM"),sep = "_"),
              paste(seq(14,42,7),c("VNIR_Plant_NDVI.avg"),sep = "_"))

h2_df = h2_df[paste(h2_df$DAT,h2_df$Trait,sep="_")%in%scenarios,]

colnames(h2_df)[2] = "Dat"
h2_df$Dat[h2_df$Trait == "RGB1_Plant_Avg_HEIGHT_MM"] = h2_df$Dat[h2_df$Trait == "RGB1_Plant_Avg_HEIGHT_MM"] - 1 

h2_df$xpos <- -Inf  # Left side
h2_df$ypos <- Inf   # Top side

plotDF = merge(plotDF,h2_df)
plotDF$CorrectedPA = plotDF$PredictionAccuracy * sqrt(plotDF$h2)
plotDF$CorrectedPA[plotDF$Model=="GBLUP"] = plotDF$PredictionAccuracy[plotDF$Model=="GBLUP"] 

plotDF=plotDF[plotDF$Trait!= "SC_Plant_Weight",]


ggplot(data = plotDF,aes(x= Model,y = CorrectedPA,fill=factor(MFE))) +
  geom_violin() +
  facet_grid(Dat~Trait) +
  geom_violin (scale="width") +
  geom_hline(yintercept = 0,linetype = "dotted")+
  geom_boxplot(notch = F,
               width = 0.1,
               position = position_dodge(width=.9)) +
  labs(fill='MBC') +
  theme_linedraw() +
  geom_text(
    data = h2_df,
    aes(x = 2, y = ypos, label = parse(text=paste("h^2:",round(h2,2)))),
    inherit.aes = F,
    hjust = -0.5,  # Nudge right
    vjust = 1.3,   # Nudge down
    color = "black",
    fontface = "bold"
  )

paste("h2:",h2_df$h2)

ggsave("../Figures/GP_UV_old.png",width = 10,height = 5)

# ---- Multivariate GP ---- 

acc_files = list.files("../Data/Generated/GenomicPrediction/MVGP/Accuracy/")
pred_files = list.files("../Data/Generated/GenomicPrediction/MVGP/Predictions/")


MV_accuracy_df = data.frame()
MV_prediction_df = data.frame()
MV_FE_df = data.frame()

acc_nf = T
pred_nf= T
fe_nf = T

for(file in acc_files){
  if(grepl("_CV2_",file)){
    CV2 = T
  }else{
    CV2 = F
  }
  if(grepl("_MFE_",file)){
    MFE = T
  }else{
    MFE = F
  }
  if(MFE){
    if(CV2){
      scenario = sub("MFE_CV2_accuracy.csv","",file)
    }else{
      scenario = sub("MFE_accuracy.csv","",file)
    }
    
  }else{
    if(CV2){
      scenario = sub("_CV2_accuracy.csv","",file)
    }else{
      scenario = sub("_accuracy.csv","",file)
    }
    
  }
  t = strsplit(scenario,"_")
  dat = t[[1]][length(t[[1]])]
  trait = paste0((t[[1]][-length(t[[1]])]),collapse = "_")
  df = read.csv(paste0("../Data/Generated/GenomicPrediction/MVGP/Accuracy/",file))
  #if(MFE){
  #  colnames(df)[c(3,2)] = c("Uhat_accuracy","Eta_mean_accuracy")
  #}else{
  #  colnames(df)[c(2,3)] = c("MegaLMM","MegaLMM_secondary")
  #}
  df$MFE = MFE
  df$CV2 = CV2
  if(acc_nf){
    MV_accuracy_df = df
    acc_nf = F
  }else{
    MV_accuracy_df = bind_rows(MV_accuracy_df, df)
  }
}
for(file in pred_files){
  if(grepl("_prediction", file)){
    if(grepl("_CV2_",file)){
      CV2 = T
    }else{
      CV2 = F
    }
    if(grepl("_MFE_",file)){
      MFE = T
    }else{
      MFE = F
    }
    if(MFE){
      if(CV2){
        scenario = sub("MFE_CV2_prediction.csv","",file)
      }else{
        scenario = sub("MFE_prediction.csv","",file)
      }
      
    }else{
      if(CV2){
        scenario = sub("_CV2_prediction.csv","",file)
      }else{
        scenario = sub("_prediction.csv","",file)
      }
      
    }
    t = strsplit(scenario,"_")
    dat = t[[1]][length(t[[1]])]
    trait = paste0((t[[1]][-length(t[[1]])]),collapse = "_")
    df = read.csv(paste0("../Data/Generated/GenomicPrediction/MVGP/Predictions/",file),row.names = 1)
    df$MFE = MFE
    df$CV2 = CV2
    if(pred_nf){
      MV_prediction_df = df
      pred_nf = F
    }else{
      MV_prediction_df = bind_rows(MV_prediction_df, df)
    }
  }
  else if(grepl("_fe_frame", file)){
    scenario = sub("_fe_frame_MegaLMM_MFE.csv","",file)
    t = strsplit(scenario,"_")
    dat = t[[1]][length(t[[1]])]
    trait = paste0((t[[1]][-length(t[[1]])]),collapse = "_")
    df = read.csv(paste0("../Data/Generated/GenomicPrediction/MVGP/Predictions/",file),row.names = 1)
    df$Trait = trait
    df$Dat = dat
    if(fe_nf){
      MV_FE_df = df
      fe_nf = F
    }else{
      MV_FE_df = bind_rows(MV_FE_df, df)
    }
  }
}

MV_acc_MFE = MV_accuracy_df %>%
  filter(MFE) %>%
  dplyr::select(DAT,Trait,rep,MFE,MegaLMM,CV2) %>%
  filter(!Trait %in% "SC_Plant_Weight")

colnames(MV_acc_MFE)[5] = "Accuracy"

MV_acc_noMFE = MV_accuracy_df %>%
  filter(!MFE) %>%
  dplyr::select(DAT,Trait,rep,MFE,U_hat_accuracy,CV2) %>%
  filter(!Trait %in% "SC_Plant_Weight") 

colnames(MV_acc_noMFE)[5] = "Accuracy"

MV_acc_merge = rbind(MV_acc_MFE,MV_acc_noMFE)

df = MV_acc_merge %>% group_by(Trait,DAT,MFE,CV2) %>% summarize(mean=mean(Accuracy),
                                                        sd = sd(Accuracy))

## Correcting accuracy of MegaLMM via cor_g(y,y_hat) * sqrt(h^2(u_hat))

K_cut = read.csv("../Data/Genotype/B1K_GRM_red.csv",row.names = 1)
trait = "RGB1_Plant_Avg_HEIGHT_MM"
dat = 21
replicate = 1
cv2 = TRUE
mfe = FALSE
data = MV_prediction_df %>% 
  filter(Trait == trait) %>%
  filter(DAT == dat) %>%
  filter(rep == replicate) %>%
  filter(CV2 == cv2) %>%
  filter(MFE == mfe) %>%
  filter(Class == "Test") %>%
  dplyr::select(c("Geno","BLUPs","MegaLMM"))


data_long = data %>% pivot_longer(cols = c("BLUPs","MegaLMM"),values_to = "Value",names_to = "Source")

fit <- mmer(
  Value ~ 1,# fixed effects
  random = ~ vsr(usr(Source), Geno, Gu = as.matrix(K_cut)),  # A = kinship/GRM matrix
  rcov   = ~ vsr(dsr(Source), units),
  data   = data_long,
  tolParInv = 30   # tolerance for matrix inversion
)
  
Va1   <- fit$sigma$`BLUPs:Geno`
Va2   <- fit$sigma$`MegaLMM:Geno`
Cov_a <- fit$sigma$`MegaLMM:BLUPs:Geno`

r_g <- Cov_a / sqrt(Va1 * Va2)
print(r_g)

E1 <- fit$sigma$`BLUPs:units`
E2 <- fit$sigma$`MegaLMM:units`

h_2_t1 = Va1 / (Va1+E1)
h_2_t2 = Va2 / (Va2+E2)



ggplot(data = MV_acc_merge,aes(x=MFE,y = Accuracy,fill=factor(CV2))) +
  geom_violin() +
  facet_grid(DAT~Trait) +
  geom_violin (scale="width")+
  geom_hline(yintercept = 0,linetype = "dotted")+
  geom_boxplot(notch = F,
               width = 0.1,
               position = position_dodge(width=.9)) +
  labs(fill='Model') +
  theme_linedraw() +
  geom_text(
    data = h2_df,
    aes(x = 2, y = -0.5, label = parse(text=paste("h^2:",round(h2,2)))),
    inherit.aes = F,
    hjust = -0.5,  # Nudge right
    vjust = 1.3,   # Nudge down
    color = "black",
    fontface = "bold"
  )



colnames(MV_accuracy_df)[6] = "Dat" 


MV_accuracy_df_n = MV_acc_merge %>%
  pivot_wider(values_from = Accuracy,names_from = CV2)

colnames(MV_accuracy_df_n)[c(1,3,5,6)] = c("Dat","Repetition","MegaLMM_CV1","MegaLMM") 

UV_accuracy_df = UV_accuracy_df %>%
  filter(Trait != "SC_Plant_Weight")

merged_acc = merge(UV_accuracy_df,MV_accuracy_df_n)

merged_acc

plotDF <- merged_acc %>%
  dplyr::select(GBLUP,HBLUP,GHBLUP,MegaLMM,MegaLMM_CV1,Trait,Dat,MFE,Repetition) %>%
  pivot_longer(cols = c(GBLUP,HBLUP,GHBLUP,MegaLMM,MegaLMM_CV1),names_to = "Model",values_to = "PredictionAccuracy")


h2_df = read.csv("../Data/Generated/nonHSR_h2.csv",row.names = 1)

scenarios = c(paste(seq(15,43,7),c("RGB1_Plant_Avg_HEIGHT_MM"),sep = "_"),
              paste(seq(14,42,7),c("VNIR_Plant_NDVI.avg"),sep = "_"))

h2_df = h2_df[paste(h2_df$DAT,h2_df$Trait,sep="_")%in%scenarios,]

colnames(h2_df)[2] = "Dat"
h2_df$Dat[h2_df$Trait == "RGB1_Plant_Avg_HEIGHT_MM"] = h2_df$Dat[h2_df$Trait == "RGB1_Plant_Avg_HEIGHT_MM"] - 1 

h2_df$xpos <- -Inf  # Left side
h2_df$ypos <- Inf   # Top side





plotDF = merge(plotDF,h2_df)
plotDF$CorrectedPA = plotDF$PredictionAccuracy * sqrt(plotDF$h2)


ggplot(data = plotDF,aes(x= MFE,y = CorrectedPA,fill=factor(Model,levels = c("GBLUP","HBLUP","GHBLUP","MegaLMM","MegaLMM_CV1")))) +
  geom_violin() +
  facet_grid(Dat~Trait) +
  geom_violin (scale="width")+
  geom_hline(yintercept = 0,linetype = "dotted")+
  geom_boxplot(notch = F,
               width = 0.1,
               position = position_dodge(width=.9)) +
  labs(fill='Model',y = "Corrected prediction accuracy",x="Marker based covariates included.") +
  theme_linedraw() +
  scale_fill_manual(values = c("#e02b35","#f0c571","#59a89c","#a559aa","#cecece"),breaks=c("GBLUP","HBLUP","GHBLUP","MegaLMM","MegaLMM_CV1")) +
  geom_text(
    data = h2_df,
    aes(x = 1, y = 0.6, label = parse(text=paste("h^2:",round(h2,2)))),
    inherit.aes = F,
    hjust = 1,  # Nudge left
    vjust = 0,   # Nudge up
    color = "black",
    fontface = "bold"
  )


ggsave("../Figures/Full_GP_accuracies.png",height=7,width=12)

write.csv(plotDF,"../Figures/GP_accuracy_data.csv")

paste("h2:",h2_df$h2)




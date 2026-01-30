## ---------------------------
##
## Script name: LD_decay_HPC.R
##
## Purpose of script: Compute and plot chromosome-wise LD decay distance, 
##                    given TASSEL-generated LD tables for each chromosome. 
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
##        The model used is Hill & Weir (1988) / Remington (2001) 
##   
##
## ---------------------------

## load up the packages we will need:  
library(ggplot2)
library(ggpubr)
library(minpack.lm)

## ---------------------------

# arg1 : Input folder
# arg2 : Sample proportion
# arg3 : Distance where no LD is expected
# arg4 : Output folder


args = commandArgs(trailingOnly = T)

plot_list=list()

outfile=paste0(args[4],"/LD_parameters.tab")

cat("Chromosome\tBG_LD\tsam_size\tHW_c",file=outfile,append=T)

for(chr in 1:7){

  
# filter (p-value < 0.05 ) and order LD table
LD_table_ch=read.table(paste0(args[1],"/B1K_LD_chr",chr,"_full.txt"),header=T)
LD_table_ch_filt = LD_table_ch[as.numeric(LD_table_ch$pDiseq)<0.05,]
LD_table_ch_filt = LD_table_ch_filt[order(LD_table_ch_filt$Dist_bp,decreasing=F),]


# calculate rolling mean (window-size 100 markers)
if(!file.exists(paste0(args[4],"/rolling_mean_chr",chr,".Rds"))){
roll_mean=c()
roll_dist=c()

for(i in 50:(dim(LD_table_ch_filt)[1]-50)){
  roll_mean=c(roll_mean,mean(na.omit(LD_table_ch_filt$R.2[(i-49):(i+50)])))
  roll_dist=c(roll_dist,LD_table_ch_filt$Dist_bp[i])
}

rolling_mean_df = data.frame(Dist=roll_dist,R2=roll_mean)
print("Rolling mean calculated")
saveRDS(rolling_mean_df,file=paste0(args[5],"rolling_mean_chr",chr,".Rds"))
}else{
  print("Rolling mean loaded")
  rolling_mean_df = readRDS(paste0(args[5],"rolling_mean_chr",chr,".Rds"))
}


# 'Custom' LD_decay functions
# Exponential decay
LD_decay <-function(t,l0,lam){return(l0*exp(-lam*t))}
# Exponential decay starting from 1 (l0 = 1)
LD_decay_fixl0 <- function(t,lam){return(1*exp(-lam*t))}
# Hill & Weir 1988
LD_decay_HW <- function(d,C){return(1/(1+C*d))}
# Remington 2001
LD_decay_HW_adj <- function(d,c,n){
  C = c*d
  r = ( (10+C) / ( (2+C) * (11+C) ) ) * ( 1 + ((3 + C) * (12 + 12*C + C^2) ) / (n * (2 + C) * (11 + C)) )
  return(r)
}

# full dataset
data_full = data.frame(Dist=as.numeric(LD_table_ch_filt$Dist_bp),R2 = as.numeric(LD_table_ch_filt$R.2))

# sample dataset
sam_size = floor(dim(data_full)[1]*as.numeric(args[2]))
data_sam = na.omit(data_full[sample(nrow(data_full),sam_size,replace=F),])

## Fit models to sample data
# Exponential decay
fit_decay         <- nlsLM(R2~LD_decay(Dist,l0,lam),data = data_sam,start = list(l0=1,lam=0.1))
# Exponential decay starting from 1 (l0 = 1)
fit_decay_fixl0   <- nlsLM(R2~LD_decay_fixl0(Dist,lam),data = data_sam,start = list(lam=0.1))
# 2nd degree loess fit
fit_decay_loess   <- loess(R2~Dist,data=data_sam,degree=2)
# Hill & Weir 1988
fit_decay_HW      <- nlsLM(R2~LD_decay_HW(Dist,C),data=data_sam,start=list(C=0.1),lower = c(0))
# Remington 2001
fit_decay_HW_adj  <- nlsLM(R2~LD_decay_HW_adj(Dist,c,sam_size),data=data_sam,start=list(c=0.1),lower=c(0))

## Create loess data frame
data_loess = data.frame(Dist=fit_decay_loess$x,R2=fit_decay_loess$fitted)
data_loess = data_loess[order(data_loess$Dist),]

## Subset the 'uncorrelated' data for background LD estimation
data_uncorr=data_full[data_full$Dist>as.numeric(args[3]),]
background_LD = quantile(data_uncorr$R2,probs = 0.95,na.rm = T)
# Create 20 quantiles of LD distribution in 'uncorrelated' pairs
background_LD_array = data.frame(BG_LD=quantile(data_uncorr$R2,probs=seq(0,1,0.05),na.rm = T))


## Compute the LD-decay at different background LD quantiles
# Exponential decay
LD_cutfreel0  <- function(x,fit) {return(log(x/coef(fit)['l0'])/(-coef(fit)['lam']))}
background_LD_array$freel0_cutoff=LD_cutfreel0(background_LD_array$BG_LD,fit_decay)
# Exponential decay starting from 1 (l0 = 1)
LD_cutfixl0   <- function(x,fit) {return(log(x/1)/(-coef(fit)['lam']))}
background_LD_array$fixl0_cutoff=LD_cutfixl0(background_LD_array$BG_LD,fit_decay_fixl0)
# Loess fit
LD_cutloess <- function(x,data_loess){return(max(data_loess$Dist[data_loess$R2>x],na.rm = T))}
background_LD_array$loess_cutoff= sapply(X=background_LD_array$BG_LD,FUN = LD_cutloess,data_loess=data_loess)
# Hill & Weir (1988)
LD_cutHW <-function(x,fit) {return((-x+1)/(x*coef(fit)['C']))}
background_LD_array$HW_cutoff = LD_cutHW(background_LD_array$BG_LD,fit_decay_HW)
# Remington (2001)
LD_cutHW_adj <- function(d,c,r,n){
  C <- c * d
  
  # First part of the equation
  part1 <- (10 + C) / ((2 + C) * (11 + C))
  
  # Second part of the equation
  part2 <- (1 + ((3 + C) * (12 + 12 * C + C^2)) / (n * (2 + C) * (11 + C)))
  
  # Full equation
  result <- part1 * part2 - r
  
  return(result)
}
c_value = coef(fit_decay_HW_adj)['c']
r_value = background_LD_array$BG_LD
n_value = sam_size
ad_cut=c()
for(i in 1:length(r_value)){
  interval_start <- 1
  interval_end <- 1e8
  
  start_value <- LD_cutHW_adj(interval_start, c = c_value, r = r_value[i], n = n_value)
  end_value <- LD_cutHW_adj(interval_end, c = c_value, r = r_value[i], n = n_value)
  
  if(start_value * end_value < 0){
    ad_cut[i]=uniroot(LD_cutHW_adj, interval = c(interval_start,interval_end),c=c_value,r=r_value[i],n=n_value)$root
  } else{
    ad_cut[i]=NA
  }
}
background_LD_array$HW_adj_cutoff = ad_cut


HW_adj_df = data.frame(Dist = seq(1,1e8,length.out=100000))
R2=c()
for(i in 1:nrow(HW_adj_df)){
  R2[i]=LD_decay_HW_adj(HW_adj_df$Dist[i],coef(fit_decay_HW_adj)['c'],sam_size)
}
HW_adj_df$R2 = R2


p1<-ggplot(data_full,aes(x=Dist,y=R2)) + 
  geom_point(size=.5,alpha=0.5) +
  geom_line(data=rolling_mean_df,aes(color="rolling mean")) +
  geom_line(data = HW_adj_df, aes(color="HW_adj"), linetype="dashed", linewidth=0.5) +
  geom_hline(yintercept = background_LD,color = "purple",linetype  ="dotted") +
  scale_color_manual("Model",values=c("green","blue","red"),breaks=c("HW_adj","rolling mean")) +
  ylim(0,1) +
  scale_x_continuous(limits = c(0,7e+08)) +
  ylab(expression("r"^2)) +
  xlab("") +
  ggtitle(paste0("Chromosome ",chr)) +
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        plot.margin = unit(c(0,0.2,0,1), 'lines'))


  
  plot_list[[chr]]=p1
  cat(paste(chr,background_LD,sam_size,coef(fit_decay_HW_adj)['c'],sep="\t"),file=outfile,append=T)
  cat("\n",file=outfile,append=T)
  write.csv(background_LD_array,file = paste0(args[4],"/LD_cutoffs_chr",chr))
}



fig=ggarrange(plot_list[[1]],
              plot_list[[2]],
              plot_list[[3]],
              plot_list[[4]],
              plot_list[[5]],
              plot_list[[6]],
              plot_list[[7]] + xlab("Dist [bp]"),
              ncol = 1,nrow=7,align="hv",common.legend = T,legend = "right") 

ggsave(paste0(args[4],"/LD_decay.png"), plot = fig, device = "png", dpi = 1200, width = 20, height = 40, units = "cm")




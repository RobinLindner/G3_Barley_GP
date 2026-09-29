## ---------------------------
##
## Script name: 
##
## Purpose of script: New LD-decay:
##        Compute LD for different MAF thresholds 
##
## Author: M.Sc. Robin Lindner
##
## Date Created: 2026-02-19
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

setwd("~/Documents/Arbeit/")    

## ---------------------------

## load up the packages we will need:  


## ---------------------------

## load up our functions into memory:

# This should be followed by the respective MAF bin upper threshold for which the LD was calculated (i.e. 0.1,0.2,0.3,0.4,and 0.5, see data directory for example)
LD_table_prefix = "../Data/Genotype/B1K_red_LD_MAF" # generated in TASSEL 

## ---------------------------


# Remington 2001
LD_decay_HW_adj <- function(d,c,n){
  C = c*d
  r = ( (10+C) / ( (2+C) * (11+C) ) ) * ( 1 + ((3 + C) * (12 + 12*C + C^2) ) / (n * (2 + C) * (11 + C)) )
  return(r)
}

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

MAF_thresholds = c(0.1,0.2,0.3,0.4,0.5)
HW_adj_curve = data.frame(Dist = seq(1,7e8,length.out=100000))
nf=T
plot_list = list()
for(MAF in MAF_thresholds){
  LD_tab = read.table(paste0(LD_table_prefix,MAF,".txt"),header = T)
  data_full = data.frame(Dist=as.numeric(LD_tab$Dist_bp),R2 = as.numeric(LD_tab$R.2)) %>%
    filter(!is.na(Dist)) %>%
    filter(!is.na(R2))
  
  sam_size = nrow(data_full)
  
  fit_decay_HW_adj  <- nlsLM(R2~LD_decay_HW_adj(Dist,c,sam_size),data=data_full,start=list(c=0.1),lower=c(0))
  
  
  data_uncorr=data_full[data_full$Dist>50e6,]
  background_LD = quantile(data_uncorr$R2,probs = 0.90,na.rm = T)
  
  # Create 20 quantiles of LD distribution in 'uncorrelated' pairs
  background_LD_array = data.frame(BG_LD=quantile(data_uncorr$R2,probs=seq(0,1,0.05),na.rm = T))
  
  # Remington (2001)
  c_value = coef(fit_decay_HW_adj)['c']
  r_value = background_LD_array$BG_LD
  n_value = sam_size
  ad_cut=c()
  for(i in 1:length(r_value)){
    interval_start <- 1
    interval_end <- 1e9
    
    start_value <- LD_cutHW_adj(interval_start, c = c_value, r = r_value[i], n = n_value)
    end_value <- LD_cutHW_adj(interval_end, c = c_value, r = r_value[i], n = n_value)
    
    if(start_value * end_value < 0){
      ad_cut[i]=uniroot(LD_cutHW_adj, interval = c(interval_start,interval_end),c=c_value,r=r_value[i],n=n_value)$root
    } else{
      ad_cut[i]=NA
    }
  }
  background_LD_array$HW_adj_cutoff = ad_cut
  
  colnames(background_LD_array) = c(paste0("BG_",MAF),paste0("CO_",MAF))
  
  cur_df = data.frame(MAF_bin = MAF,
                      LD_decay = background_LD_array["90%",2],
                      BG_LD = background_LD_array["90%",1],
                      estimated_c = c_value,
                      n = sam_size)
  if(nf){
    full_BG_LD = background_LD_array
    res_df = cur_df
    nf=F
  }else{
    full_BG_LD = cbind(full_BG_LD,background_LD_array)
    res_df = rbind(res_df,cur_df)
  }
  
  
  ## Plotting
  R2=c()
  for(i in 1:nrow(HW_adj_curve)){
    R2[i]=LD_decay_HW_adj(HW_adj_curve$Dist[i],coef(fit_decay_HW_adj)['c'],sam_size)
  }
  HW_adj_curve = HW_adj_curve %>% 
    mutate(!!toString(MAF) := R2)
  
  plotdf = data.frame(Dist = HW_adj_curve$Dist,
                      R2 = R2)
  
  p1<-ggplot(data_full,aes(x=Dist,y=R2)) + 
    geom_point(size=.5,alpha=0.5) +
    geom_line(data = plotdf, aes(color="HW_adj"), linetype="dashed", linewidth=0.5) +
    geom_hline(yintercept = background_LD,color = "purple",linetype  ="dotted") +
    geom_vline(xintercept =ad_cut[19],color = "darkorange") +
    scale_color_manual("Model",values=c("green"),breaks=c("HW_adj")) +
    ylim(0,1) +
    geom_text(x=ad_cut[19] + 5e7,y = 0.05,color = "darkorange",label = paste0(toString(round(ad_cut[19])),"bp"),inherit.aes = F)+
    scale_x_continuous(limits = c(0,7e+08)) +
    ylab(expression("r"^2)) +
    xlab("") +
    ggtitle(paste0("MAF: [",MAF-0.1,"-",MAF,")")) +
    theme_linedraw()+
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8),
          plot.margin = unit(c(0,0.2,0,1), 'lines'))
    
    
  
  
  
  plot_list[[toString(MAF)]]=p1
  
}


fig=ggarrange(plot_list[["0.1"]],
              plot_list[["0.2"]],
              plot_list[["0.3"]],
              plot_list[["0.4"]],
              plot_list[["0.5"]] + xlab("Dist [bp]"),
              ncol = 1,nrow=5,align="hv",common.legend = T,legend = "right") 

ggsave("../Supplements/LD_decay_figure_new.png", plot = fig, device = "png", dpi = 1200, width = 20, height = 40, units = "cm")



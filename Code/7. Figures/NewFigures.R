#New Figures
source("0_utils.R")
library(scales)

#### Data prep ####

remap <- read.csv(geno_remap_file,row.names = 1)

groups = read.csv(trait_groups_file)
group_levels = c("CC","HL.LL","HL","CL","VNIR","IR","RGB","SW","H")


#### I. Kinship #### 
K=read.csv("../Data/Genotype/B1K_IPK/B1K_IPK_final.VRKin.csv",row.names = 1)

full=read.csv("../Data/Phenotype/Capitalize B1K data_Tier 1-Tier 2_PSI_updated/Merged_file_Tier1_Enviro.csv")

geno_loc <- full %>%
  dplyr::select(Genotype,Location) %>%
  unique()

pca = prcomp(K)

df=data.frame(Genotype=rownames(K),Location=geno_loc$Location[match(rownames(K),geno_loc$Genotype)],PC1=pca$rotation[,1],PC2=pca$rotation[,2],PC3=pca$rotation[,3])



# Extract variance explained by each PC
variance_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
cumulative_variance <- cumsum(variance_explained)

# Number of PCs to display in the scree plot
num_pcs <- 10  # Adjust this to the desired number of PCs
num_pcs <- min(num_pcs, length(variance_explained))  # Ensure it doesn't exceed total PCs

# Create a data frame for plotting
scree_df <- data.frame(
  PC = 1:num_pcs,
  Variance_Explained = variance_explained[1:num_pcs],
  Cumulative_Variance = cumulative_variance[1:num_pcs]
)

p1<-ggplot(df,aes(x=PC1,y=PC2,color=Location))+
  geom_point(size=1) +
  xlab(paste0("PC1 [",round(variance_explained[1],digits = 2),"%]")) +
  ylab(paste0("PC2 [",round(variance_explained[2],digits = 2),"%]")) +
  scale_colour_discrete(name="") +
  labs(tag="B") +
  theme_linedraw()+
  theme(plot.tag = element_text())

p2<-ggplot(df,aes(x=PC1,y=PC3,color=Location))+
  geom_point(size=1) +
  xlab(paste0("PC1 [",round(variance_explained[1],digits = 2),"%]")) +
  ylab(paste0("PC3 [",round(variance_explained[3],digits = 2),"%]")) +
  scale_colour_discrete(name="")+
  labs(tag="C") +
  theme_linedraw()+
  theme(plot.tag = element_text())



scree_df$Included = "No"
scree_df$Included[scree_df$Variance_Explained>5] = "Yes"




p3 <-ggplot(scree_df, aes(x = PC, y = Variance_Explained,fill=factor(Included,levels=c("Yes","No")))) +
  geom_bar(stat = "identity", alpha = 1) +
  scale_x_continuous(breaks = 1:num_pcs) +
  labs(
    title = "",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  labs(tag="A") +
  scale_fill_manual(values = c("Yes"="lightgreen","No"="grey"),name="Considered in GWAS")+
  #geom_vline(xintercept = 4.5,linetype = "dotted",color = "black")+
  theme_linedraw()+
  theme(plot.tag = element_text(),
        legend.position = "bottom")

col2=ggarrange(p1,p2,nrow=2,ncol=1, common.legend = TRUE, legend="right")
ggarrange(p3,col2,nrow=1,ncol=2)

ggsave("../Figures/NewKinship_Figure.png",device = "png",bg = "white",width=12,height=5)

#### II. Heritability ####



h2_long_df=read.csv(nonHSR_h2_path,row.names = 1)
h2_long_df$Group=groups$Group[match(h2_long_df$Trait,groups$Trait)]

h2_long_df$Group[h2_long_df$Group=="TCI"] = "IR"
h2_long_df$Group = factor(h2_long_df$Group,levels = group_levels)

## Medians
medians = h2_long_df %>% 
  group_by(Group,DAT) %>%
  summarize(median=median(h2),
            sd=sd(h2),
            quant25=quantile(h2,probs = 0.25),
            quant75=quantile(h2,probs = 0.75))

#fit regression slope
counts = medians %>%
  group_by(DAT) %>%
  summarise(group_c=n())
sel_dats=counts$DAT[counts$group_c>=4]

medians_f <- medians%>%
  filter(DAT %in% sel_dats)
reg = lm.fit(x=as.matrix(cbind(rep(1,length(medians_f$DAT)),medians_f$DAT)),y=medians_f$median)

tr_n= c(6, 9, 41, 41, 6, 2, 11, 2, 3)
tp_n= c(6, 19, 9, 9, 6, 19, 35, 40, 1)



label_size = 5
text_size = 14

p1<-ggplot(h2_long_df,aes(x=Group,y=h2,fill=Group)) + 
  geom_violin(scale = "width") +
  geom_boxplot(notch = F,
               width = 0.1) +
  labs(y=expression(h^2),
       x="Measurement protocol",
       tag="A") +
  # annotate("label",
  #          x = seq(1,9),
  #          y = 0.78,
  #          label = tr_n,
  #          hjust = 0.5,
  #          vjust = -1.7,
  #          size = label_size,
  #          color = "black",
  #          fill="lightgrey") +
  # annotate("label",
  #          x = seq(1,9),
  #          y = 0.78,
  #          label = tp_n,
  #          hjust = 0.5,
  #          vjust = -0.5,
  #          size = label_size,
  #          color = "black",
  #          fill="lightblue") +
  scale_fill_discrete(breaks =levels(h2_long_df$Group),labels= trait_groups) +
  lims(y=c(0,0.8))+
  theme_linedraw()+
  theme(legend.position="none",
        axis.text.x = element_text(size = text_size),
        axis.text.y = element_text(size = text_size),
        axis.title = element_text(size = text_size),
        legend.title = element_text(size = text_size*1.2),
        legend.text = element_text(family = "mono",face = "bold",size = text_size))

ggsave("../Figures/Heritability_Figure_plot1.png",plot=p1,height = 6,width=14,device = "png")

trait_groups = c(
  "CC:    Chlorophyll content",
  "HL.LL: High-light/low-light shift response",
  "HL:    High light PSII efficiency",
  "CL:    Cultivation light PSII efficiency",
  "VNIR:  Vitality Indices",
  "IR:    Thermoregulation",
  "RGB:   Morphological features",
  "SW:    Weight traits",
  "H:     Biomass at Harvest"
)
trait_slopes =c(
  "Chlorophyll content:                 −6.3e−3",
  "High-light/low-light shift response: +2.5e-3",
  "High light PSII efficiency:          +8.3e−3",
  "Cultivation light PSII efficiency:   +6.1e−3",
  "Vitality Indices:                    +0.014",
  "Thermoregulation:                    +2.0e-3",
  "Morphological features:              +3.3e-3",
  "Weight traits:                       +4.3e-4",
  "Biomass at Harvest:                   ------"
)


medians$linetype = "dotted"
medians$linetype[medians$Group %in%c("RGB","VNIR")] = "solid"
p2<-ggplot(medians,aes(x=DAT,y=median,color=Group)) +
  geom_abline(color="darkgray",
              slope = reg$coefficients[2],
              intercept = reg$coefficients[1],
              linetype="longdash",
              linewidth=1)+
  #geom_errorbar(aes(ymin=median-sd,ymax=median+sd),width=.1,linetype = "dotted") +
  geom_line(linetype = medians$linetype,linewidth = .7, show.legend = F) +
  geom_point() +
  scale_color_discrete(name="Average change in heritability per day",breaks =levels(medians$Group),labels= trait_slopes) +
  labs(y=expression("Median"~h^2),
       x="Days after transplant",
       tag="B") +
  lims(y=c(0,0.8))+
  theme_linedraw()+
  theme(legend.position="none",
    axis.text.x = element_text(size = text_size),
    axis.text.y = element_text(size = text_size),
    axis.title = element_text(size = text_size),
    legend.title = element_text(size = text_size*1.2),
    legend.text = element_text(family = "mono",face = "bold",size = text_size))

ggarrange(p1,p2,ncol = 2,align = "h")
ggsave("../Figures/Heritability_Figure_newest.png",height=5,width=16,device = "png")


#### III. Pleiotropic markers ####

all_assoc = read.csv(sig_associations_file)

snp_48_traits = "chr4H_632274504"

snp_37_traits = "chr2H_704876309"

snp_30_traits = "chr1H_349277204"

#snp_29_traits = "chr5H_621066587"

### find new positions
new_pos <- remap %>%
  filter(SNP %in% c(snp_48_traits,snp_37_traits,snp_30_traits))

### Compare trait sets of the markers

traits1 <- all_assoc %>%
  filter(SNP==snp_48_traits) %>%
  dplyr::select(DAT,Trait) %>%
  distinct(.keep_all = T) %>%
  mutate(Group = trait_groups$Group[match(Trait,trait_groups$Trait)]) %>%
  group_by(Trait,Group) %>%
  summarize(nr_tp=n())

length(unique(traits1$Group))

traits2 <- all_assoc %>%
  filter(SNP==snp_37_traits) %>%
  dplyr::select(DAT,Trait) %>%
  distinct(.keep_all = T) %>%
  mutate(Group = trait_groups$Group[match(Trait,trait_groups$Trait)]) %>%
  group_by(Trait,Group) %>%
  summarize(nr_tp=n())

length(unique(traits2$Group))


traits3 <- all_assoc %>%
  filter(SNP==snp_30_traits) %>%
  dplyr::select(DAT,Trait) %>%
  distinct(.keep_all = T) %>%
  mutate(Group = trait_groups$Group[match(Trait,trait_groups$Trait)]) %>%
  group_by(Trait,Group) %>%
  summarize(nr_tp=n())

length(unique(traits3$Group))

## Panel A - all associations ##

assoc_cut_noDAT <- all_assoc %>%
  select(SNP,Trait,DAT) %>%
  distinct(.keep_all = T) %>%
  mutate(Group = trait_groups$Group[match(Trait,trait_groups$Trait)]) %>%
  select(SNP,Trait,Group)%>%
  distinct(.keep_all=T)


snp_trait_counts <- assoc_cut_noDAT %>%
  group_by(SNP) %>%
  summarise(Count=n()) %>%
  mutate(Bin = cut(Count,breaks=seq(1,max(Count),length.out=),include.lowest = T,labels = F))

assoc_cut_noDAT_wbin <- assoc_cut_noDAT %>%
  mutate(trait_n = snp_trait_counts$Count[match(assoc_cut_noDAT$SNP,snp_trait_counts$SNP)])

sum_per_trait_n = assoc_cut_noDAT_wbin %>%
  group_by(trait_n) %>%
  summarise(Count=length(unique(SNP)))

plot_df <- assoc_cut_noDAT_wbin %>%
  group_by(Group,trait_n) %>%
  summarize(Count=length(unique(SNP)))


plot_df$Group[plot_df$Group=="TCI"] = "IR"

plot_df_sums <- plot_df %>%
  group_by(trait_n) %>%
  summarise(sum = sum(Count))

plot_df <- plot_df %>%
  mutate(sum = plot_df_sums$sum[match(trait_n,plot_df_sums$trait_n)],
         snps = sum_per_trait_n$Count[match(trait_n,sum_per_trait_n$trait_n)]) %>%
  mutate(prop=Count/sum,
         adj_snps = (Count/sum)*snps)


snp_trait_nrs = all_assoc %>%
  group_by(SNP) %>%
  summarise(uTraits = length(unique(Trait)))



group_levels = c("CC","HL.LL","HL","CL","VNIR","IR","RGB","SW","H")
trait_group_names =c(
  "Chlorophyll content",
  "High-light/low-light shift response",
  "High light PSII efficiency",
  "Cultivation light PSII efficiency",
  "Vitality Indices",
  "Thermoregulation",
  "Morphological features",
  "Weight traits",
  "Biomass at Harvest"
)
group_colors = setNames(hue_pal()(length(group_levels)), group_levels)

text_size = 12

p1<-ggplot(plot_df,aes(x=trait_n,y=adj_snps,fill=factor(Group,levels=group_levels))) +
  geom_bar(position = "stack",stat="identity") +
  geom_text(data=sum_per_trait_n,
             aes(x=trait_n,y=Count,label = Count),
             vjust = -0.5, 
             size = text_size/3,
             inherit.aes = F) + # Add labels above the bars
  geom_label(data=data.frame(x_d=c(30,37,48),y_d=350,label_d=c("pSNP3","pSNP2","pSNP1")),
             aes(x=x_d,y=y_d,label=label_d),
             angle = 90,
             inherit.aes = F)+
  scale_fill_manual(values=group_colors,labels = trait_group_names) +
  scale_x_continuous(
    breaks = c(1,seq(10,50,10)),
    minor_breaks = seq(1,50,1),
    guide = guide_axis(minor.ticks = TRUE)) +
  lims(y=c(0,1500))+
  labs(x="Number of associated traits",
       y="Number of markers",
       fill="Group",
       tag="A") +
  theme_linedraw()+
  theme(legend.position = c(.8,.7),
        legend.title = element_text(size = text_size+4), 
        legend.text = element_text(size = text_size,lineheight = 3),
        legend.key.size = unit(4,'mm'),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        axis.text.x = element_text(size=text_size),
        axis.title = element_text(size=text_size),
        axis.text.y = element_text(size=text_size)) 


## Panel B - top associations ##

trait_sum1 <- traits1 %>%
  group_by(Group) %>%
  summarize(Ntrait = n()) %>%
  mutate(SNP = "pSNP1")

trait_sum2 <- traits2 %>%
  group_by(Group) %>%
  summarize(Ntrait = n()) %>%
  mutate(SNP = "pSNP2")

trait_sum3 <- traits3 %>%
  group_by(Group) %>%
  summarize(Ntrait = n()) %>%
  mutate(SNP = "pSNP3")

# trait_sum4 <- all_assoc %>%
#   filter(SNP=="chr1H_186428634") %>%
#   dplyr::select(DAT,Trait) %>%
#   distinct(.keep_all = T) %>%
#   mutate(Group = trait_groups$Group[match(Trait,trait_groups$Trait)]) %>%
#   group_by(Trait,Group) %>%
#   summarize(nr_tp=n()) %>%
#   group_by(Group) %>%
#   summarize(Ntrait = n()) %>%
#   mutate(SNP = "pSNP4")
# 
# trait_sum5 <- all_assoc %>%
#   filter(SNP=="chr1H_349277204") %>%
#   dplyr::select(DAT,Trait) %>%
#   distinct(.keep_all = T) %>%
#   mutate(Group = trait_groups$Group[match(Trait,trait_groups$Trait)]) %>%
#   group_by(Trait,Group) %>%
#   summarize(nr_tp=n()) %>%
#   group_by(Group) %>%
#   summarize(Ntrait = n()) %>%
#   mutate(SNP = "pSNP5")



plot_df2 = rbind(trait_sum1,
                 trait_sum2,
                 trait_sum3)#,
                 #trait_sum4,
                 #trait_sum5)

p2<-ggplot(plot_df2,aes(x=SNP,y=Ntrait,fill = factor(Group,levels=group_levels)))+
  geom_bar(position="stack",stat="identity",width = 0.5 )+
  geom_label(
    data=data.frame(tot = c(48,37,30),pos=c(1,2,3),lab=c(48,37,30)),
    aes(x=pos,y=tot,label = lab),
    inherit.aes = F,hjust = - 0.3
  )+
  coord_flip() +
  scale_fill_manual(values=group_colors,labels = trait_group_names)+
  labs(x="",
       y="Number of associated traits",
       fill = "Group",
       tag="B")+
  lims(y=c(0,50))+
  theme_linedraw()+
  theme(legend.position = "none",
        legend.title = element_text(size = text_size+4), 
        legend.text = element_text(size = text_size,lineheight = 3),
        
        legend.key.size = unit(4,'mm'),
        axis.text.x = element_text(size=text_size),
        axis.title = element_text(size=text_size),
        axis.text.y = element_text(size=text_size)) 
  
ggarrange(p1,p2,nrow=1,ncol=2, align="h",widths = c(2,1))

ggsave("../Figures/PleiotropicMarkers_newest.png",device="png",height=5,width=16)


#### IV. PH GWAS ####

ph_map = read.csv(ph_snp_ID_map)
heightPvalue_dir = "../Data/Generated/GWAS_results/Height/"
files = list.files(heightPvalue_dir)
nf=T
for(file in files){
  dat = sub("\\.FarmCPU.csv","",sub(".+_\\._\\._","",file))
  p_val_frame = read.csv(paste0(heightPvalue_dir,file))
  p_val_frame$DAT = dat
  colnames(p_val_frame)[9] = "p_value"
  if(nf){
    full_p_frame = p_val_frame
    nf=F
  }else{
    full_p_frame = rbind(full_p_frame,
                         p_val_frame)
  }
}

map = AddAdjustedPostionToMap(read.table("../Data/Genotype/B1K_red.geno.map",header=T))

plot_df = full_p_frame %>%
  group_by(SNP,CHROM) %>%
  summarise(MIN = min(p_value,na.rm = T),
            MAX = max(p_value,na.rm = T),
            MED = median(p_value,na.rm=T)) %>%
  mutate(adj_pos = map$adjPOS[match(SNP,map$SNP)])



text_size = 12
anno_size = 4.5

quant_p=quantile(full_p_frame$p_value,probs = seq(0,1,0.01),na.rm=T)

p2 <-ggplot(plot_df) +
  geom_segment(aes(x=adj_pos,
                   xend = adj_pos,
                   y=-log10(MED),
                   yend=-log10(MIN),
                   color=as.factor(CHROM)),
               linewidth = .1)+
  #geom_point(data = long,aes(x=adj_pos,y=-log10(MED),color=CHROM),) +
  geom_hline(yintercept = -log10(0.05/1322),color="red",linetype = "dashed")+
  #geom_hline(yintercept = -log10(quant_p[2]),color="blue",linetype = "dashed")+
  #scale_color_manual(values=c("1"="blue","2"="green","3"="orange","4"="cyan","5"="purple","6"="pink","7"="yellow"))+
  scale_color_brewer(palette = "Set1") +
  labs(tag="A")+
  theme(legend.position = "none",
        panel.background = element_rect("white","black"),
        panel.grid.major = element_line("gray",linewidth = 0.1),
        panel.grid.minor=element_line("gray",linewidth = 0.05),
        plot.tag=element_text(),
        axis.text = element_text(size=text_size),
        axis.title = element_text(size=text_size*1.1)) +
  labs(x="Genome position [bp]",y="-log10(p)")+ 
  annotate("text",
           x=max(long$adj_pos),
           y= -log10(0.05/1322),
           label = paste("p =",signif(0.05/1322,3)),
           hjust = 1.01,
           vjust = 1.2,
           size = anno_size,
           color = "red")

#quant_p=quantile(long$Value,probs = seq(0,1,0.01),na.rm=T)

th1 = 0.05/1322
th2 = 0.01

plot_df$highlight = as.numeric(plot_df$MED<=th2) + as.numeric(plot_df$MED<=th1)

significant_snps <- plot_df %>%
  filter(highlight>0) %>%
  distinct(.keep_all=T)


significant_snps_1 <- plot_df %>%
  filter(highlight==1) %>%
  dplyr::select(SNP,adj_pos,MIN,MED,highlight) %>%
  distinct(.keep_all=T)

significant_snps_2<- plot_df %>%
  filter(highlight==2) %>%
  dplyr::select(SNP,adj_pos,MIN,MED,highlight) %>%
  distinct(.keep_all=T)



p3 <- ggplot(plot_df) +
  geom_segment(aes(x=adj_pos,
                   xend = adj_pos,
                   y=-log10(MED),
                   yend=-log10(MIN),
                   color=as.factor(highlight)),
               
               #alpha=(1/34),
               linewidth = .1)+
  geom_segment(data=significant_snps,
               aes(x=adj_pos,
                   xend = adj_pos,
                   y=-log10(MED),
                   yend=-log10(MIN),
                   color=as.factor(highlight),
                   alpha=1),
               linewidth = .4)+
  #geom_point(data = long,aes(x=adj_pos,y=-log10(MED),color=CHROM),) +
  geom_hline(yintercept = -log10(th1),color="red",linetype = "dashed")+
  geom_hline(yintercept = -log10(th2),color="blue",linetype = "dashed")+
  #scale_color_manual(values=c("1"="blue","2"="green","3"="orange","4"="cyan","5"="purple","6"="pink","7"="yellow"))+
  scale_color_manual(values=c("grey","blue","red")) +
  theme(legend.position = "none",
        panel.background = element_rect("white","black"),
        panel.grid.major = element_line("gray",linewidth = 0.1),
        panel.grid.minor=element_line("gray",linewidth = 0.05),
        plot.tag=element_text(),
        axis.text = element_text(size=text_size*1.1),
        axis.title = element_text(size = text_size*1.1)) +
  labs(x="Genome position [bp]",y="-log10(p)",tag="B")+

  annotate("text",
               x = significant_snps_1$adj_pos[1],# Filter for significant SNPs
               y = -log10(significant_snps_1$MIN[1]), # Adjust y to place below axis
               label = ph_map$ID[match(significant_snps_1$SNP[1],ph_map$SNP)],     
               hjust = 1.01,           # Adjust horizontal alignment
               size = anno_size,
               color = "blue",
               
) + annotate("text",
               x = significant_snps_1$adj_pos[2],# Filter for significant SNPs
               y = -log10(significant_snps_1$MIN[2]), # Adjust y to place below axis
               label = ph_map$ID[match(significant_snps_1$SNP[2],ph_map$SNP)],     
               hjust = -0.05, # Adjust horizontal alignment
               size = anno_size,
               color = "blue",
               
) + annotate("text",
                x = significant_snps_2$adj_pos,# Filter for significant SNPs
                y = -log10(significant_snps_2$MIN), # Adjust y to place below axis
                label = ph_map$ID[match(significant_snps_2$SNP,ph_map$SNP)],     
                hjust = -0.05,           # Adjust horizontal alignment
                size = anno_size,
                color = "red"
) + annotate("text",
                 x=max(long$adj_pos),
                 y= c(-log10(th1),-log10(th2)),
                 label = paste("p =",format(signif(c(th1,th2),3),scientific=T)),
                 hjust = 1.01,
                 vjust = 1.2,
                 size = anno_size,
                 color = c("red","blue"))


p3
occ_plot=read.csv("../../../Figures/Supplements/Occurence_time_height.csv")

## We are only plotting SNP that show significant association in 3 or more consecutive time points.

has_three_consecutive <- function(x) {
  if (length(x) < 3) return(FALSE)
  # Calculate differences between adjacent elements
  # A consecutive sequence has a difference of 1
  diffs <- diff(x)
  # Check if two adjacent differences are both equal to 1
  return(any(diffs[-length(diffs)] == 1 & diffs[-1] == 1, na.rm = TRUE))
}


ph_relevant_snp = full_p_frame %>%
  filter(p_value <= 0.05/1322) %>%
  mutate(DAT = as.numeric(DAT)) %>%
  arrange(DAT) %>%
  group_by(SNP) %>%
  filter(has_three_consecutive(DAT))
  
  


alt_all_p = read.csv("../Data/Generated/RegularTraits/T1/GWA_res_4PC/RGB1_Plant_Avg_HEIGHT_MM_all_p_values.csv")
alt_all_p = alt_all_p[,c(1,2,3,order(as.numeric(sub("X","",colnames(alt_all_p)[4:38])))+3)]
na_frame = which(is.na(alt_all_p),arr.ind = T)
for(i in 1:nrow(na_frame)){
  row=na_frame[i,1]
  col=na_frame[i,2]
  alt_all_p[row,col] = mean(c(alt_all_p[row,col-1],alt_all_p[row,col+1]),na.rm=T)
}

long2 = alt_all_p %>%
  pivot_longer(cols=colnames(alt_all_p)[4:38],values_to = "p_value",names_to = "DAT") %>%
  mutate(DAT = as.numeric(sub("X","",DAT)))

p_values=c()
for(i in 1:nrow(occ_plot)){
  t= long2 %>%
    filter(SNP==occ_plot$SNP[i],DAT==occ_plot$DAT[i])
  p_values = c(p_values,t$p_value)
}
occ_plot$p_value=p_values

ph_map = data.frame(SNP=unique(c(ph_relevant_snp$SNP,significant_snps$SNP)))
ph_map$chrom = substring(sub("chr","",ph_map$SNP),1,1)
ph_map$pos = sub(".+_","",ph_map$SNP)
ph_map = ph_map[with(ph_map,order(chrom,as.numeric(pos))),]

ph_map$num_on_chrom = c(seq(1,sum(ph_map$chrom==1)),
                        seq(1,sum(ph_map$chrom==2)),
                        seq(1,sum(ph_map$chrom==3)),
                        seq(1,sum(ph_map$chrom==4)),
                        seq(1,sum(ph_map$chrom==5)),
                        seq(1,sum(ph_map$chrom==7)))

ph_map$ID = paste0("ph-",ph_map$chrom,"-",ph_map$num_on_chrom) 


occ_plot = full_p_frame %>%
  filter(SNP %in% c(ph_map$SNP))

occ_plot$p_value[is.na(occ_plot$p_value)] =1

occ_plot$SNPID = ph_map$ID[match(occ_plot$SNP,ph_map$SNP)]

write.csv(ph_map,ph_snp_ID_map,row.names = F)


p4 <- ggplot(occ_plot,
             aes(factor(DAT,levels=sort(unique(as.numeric(DAT)))),
                 factor(SNPID,levels =ph_map$ID),
                 fill = -log10(p_value))) +
  geom_tile(color="black",linewidth = .2) + 
  geom_point(data = occ_plot[(occ_plot$p_value<=th2 & occ_plot$p_value>th1),],shape = 8,color="blue",alpha=0.5) +
  geom_point(data = occ_plot[occ_plot$p_value<=th1,],shape = 8,color="red",alpha=0.5) +
  theme(axis.text = element_text(size =text_size)) +
  scale_fill_gradient(low="white",high="forestgreen") +
  geom_vline(xintercept = c(13.5,24.5),colour = "black",linetype = "dashed",linewidth = 1.1) +
  annotate("text", 
           x = c(6.25,18,29), 
           y = 32,  # Slightly above plot
           vjust=-1.5,
           label = c("early","intermediate","late"), 
           size = 5,
           fontface = "bold") +
  coord_cartesian(clip = "off") +
  
  labs(x = "Days after transplant",
       y = "SNP",
       tag="C") +
  theme_linedraw() +
  theme(legend.position="none",
        plot.tag = element_text(),
        plot.margin = margin(t = 30, r = 10, b = 10, l = 10),
        axis.text.y = element_text(size=text_size*0.9),
        axis.text.x = element_text(size=text_size*0.9,
                                   angle = 45,
                                   hjust = 1,
                                   vjust = 1),
        axis.title = element_text(size=text_size*1.1))

p4 

p1 <- ggarrange(p2,p3,nrow=2,ncol=1,align="v")
p_full <- ggarrange(p1,p4,ncol=2,nrow=1,widths = c(1,1)) 

ggsave("../Figures/height_association_new.png",plot = p_full,height = 7,width = 16,device="png")


#### V. Genomic Prediction ####

all_acc_frame = read.csv("../Data/Generated/GenomicPrediction/all_trait_accuracy.csv")
all_MFE_acc_frame = read.csv("../Data/Generated/GenomicPrediction/all_trait_accuracy_MFE.csv")
comb_all_trait_acc_df=bind_rows(all_acc_frame,all_MFE_acc_frame)
comb_all_trait_acc_df$Method = c(rep("no_MFE",1500),rep("MFE",1500))


cat_long = read.csv("../Figures/GP_accuracy_data.csv",row.names = 1)
cat_long_cut = cat_long %>%
  dplyr::select(Trait,MFE,Model,CorrectedPA)
cat_long_cut$Method =NA
cat_long_cut$Method[cat_long_cut$MFE==T] = "Included"
cat_long_cut$Method[cat_long_cut$MFE==F] = "Excluded"
cat_long_cut$Model[cat_long_cut$Model=="GHBLUP"] = "G+HBLUP"
cat_long_cut$Model[cat_long_cut$Model=="MegaLMM_CV1"] = "MegaLMM CV1"
#cat_long_cut$Model[cat_long_cut$Model=="GBLUP"] = "GBLUP"
#cat_long_cut$Model[cat_long_cut$Model=="HBLUP"] = "HBLUP"

cat_long_cut$Trait[cat_long_cut$Trait=="RGB1_Plant_Avg_HEIGHT_MM"] = "Plant height [h<sup>2</sup> ≈ 0.61]"
cat_long_cut$Trait[cat_long_cut$Trait=="VNIR_Plant_NDVI.avg"] = "NDVI [h<sup>2</sup> ≈ 0.53]"
#cat_long_cut$Trait[cat_long_cut$Trait=="SC_Plant_Weight"] = "Pot weight [H<sup>2</sup> ≈ 0.02]"


dodge_width=0.9
ggplot(cat_long_cut,
       aes(x= factor(Model,levels=c("GBLUP","HBLUP","G+HBLUP","MegaLMM","MegaLMM CV1")),
           y = CorrectedPA,
           fill = factor(Method,
                         levels=c("Excluded",
                                  "Included"))))+
  geom_violin (scale="width")+
  geom_boxplot(notch = F,
               width = 0.1,
               position = position_dodge(width=.9)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 22, size = 5, color = NA))) +
  # geom_bar(stat="identity",position=position_dodge(width = dodge_width),width=dodge_width) +
  # geom_errorbar(aes(ymax=upper_bound,
  #                   ymin=lower_bound),
  #               position=position_dodge(width = dodge_width),
  #               width = 0.4) +
  
  #geom_hline(yintercept = 0,linetype = "dotted")+
  facet_wrap(~ factor(Trait,levels=c("Plant height [h<sup>2</sup> ≈ 0.61]","NDVI [h<sup>2</sup> ≈ 0.53]")),nrow=2) + 
  geom_hline(yintercept = 0,linetype="dotted")+
  scale_fill_discrete(name="MBC") +
  labs( x="", 
        y="Prediction accuracy") +
  theme_linedraw() +
  theme(legend.position = "right",
        legend.text = element_text(size=text_size),
        legend.title = element_text(size=text_size*1.1),
        strip.background = element_rect(fill="lightgrey"),
        strip.text = element_markdown(color="black",
                                  size=text_size,),
        axis.text = element_text(size=text_size),
        axis.title = element_text(size=text_size*1.1),
        #axis.text.x = element_text(angle=4,hjust=1)
        )


ggsave("../Figures/newGP_figure.png",height=7,width=12,device = "png")


### ---- Table3: PH markers ----

ph_map = read.csv(ph_snp_ID_map)
remap = read.csv(geno_remap_file)

# ID | Chromosome | Position | Developmental period | Candidate genes
newPos = remap$new_position[match(ph_map$SNP,remap$SNP)]
newPos[is.na(newPos)]=ph_map$pos[is.na(newPos)]

sites = read.table("../Data/Genotype/B1K_red_siteSummary.txt",sep="\t",header=T)
sites = sites %>%
  select(Site.Name,Chromosome,Physical.Position,Minor.Allele.Frequency)

maf_vec = sites$Minor.Allele.Frequency[match(ph_map$SNP,sites$Site.Name)]

# variable is calculated in 5_GO_enrichment
cand_gene_num = unlist(lapply(FUN=length,X = cand_genes))
cand_gene_snp = names(cand_genes)

cand_gene_vec = c()
cand_gene_vec[match(cand_gene_snp,ph_map$ID)] = cand_gene_num

early_snp = ph_map$SNP[c(1:11,18,25,28,31)]
inter_snp = ph_map$SNP[c(12,15,18,20,26,27,29,32)]
late_snp = ph_map$SNP[c(13,14,16:24,30,32)]

early_idx = c(1:11,18,25,28,31)
inter_idx = c(12,15,18,20,26,27,29,32)
late_idx = c(13,14,16:24,30,32)
e_i_inter = intersect(early_idx,inter_idx)
e_l_inter = intersect(early_idx,late_idx)
i_l_inter = intersect(inter_idx,late_idx)

period_vec = c()
period_vec[early_idx] = "Early"
period_vec[inter_idx] = "Intermediate"
period_vec[late_idx] = "Late"
period_vec[i_l_inter] = "Intermediate - Late"
period_vec[18] = "All"

ph_marker_table = data.frame(ID = ph_map$ID,Chromosome=paste0(ph_map$chrom,"H"),Position=newPos,MinorAlleleFrequency = maf_vec,DevelopmentalPeriod=period_vec,CandidateGenes=cand_gene_vec)

write.csv(ph_marker_table,"../Figures/Table3.csv")




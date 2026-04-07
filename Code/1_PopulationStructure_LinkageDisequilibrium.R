## ---------------------------
##
## Script name: 1_PopulationStructure_LinkageDisequilibrium
##
## Purpose of script: Contains analyses conducted for named aspects of the G3
##                    Barley time-series GWAS and GP publication.
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


source("0_utils.R")

## ---------------------------

## ---- Population structure ----

# Generate the rMVP genotype and genome map files
MVP.Data(fileVCF = vcf_file,out = rMVP_out)


# Compute the VanRaden genomic relatedness matrix (GRM)
genotype = attach.big.matrix(geno_path)
dim(genotype)


if(!file.exists(GRM_path)){
  K = MVP.K.VanRaden(genotype,mrk_bycol = T)
  geno_IDs = read.table(geno_IDs_path)
  rownames(K)=geno_IDs$V1
  write.csv(K,GRM_path,row.names = T)
}else{
  K = read.csv(GRM_path,row.names = 1)
}

renameVCFTaxa <- function(infile,outfile,prelim_col = 9){
  lines = readLines(infile)
  i==1
  while(!startsWith(lines[i],"#CHROM")){
    i=i+1
  }
  ID_line = lines[i]
  split = strsplit(ID_line,"\t")
  taxa_idxs = (prelim_col+1):length(split[[1]])
  oldtaxa = split[[1]][taxa_idxs]
  newtaxa = transformTaxaToGenotypeID(oldtaxa)
  split[[1]][taxa_idxs] = newtaxa
  lines[i]=paste(split[[1]],collapse="\t")
  writeLines(lines,outfile)
}

renameVCFTaxa("../Data/Genotype/B1K_red_ot.vcf","../Data/Genotype/B1K_red.vcf")

geno = as.matrix(attach.big.matrix("../Data/Genotype/B1K_red.geno.desc"))
sites = read.table("../Data/Genotype/B1K_red.geno.map",header=T)
tid = read.table("../Data/Genotype/B1K_red.geno.ind")


## Create numeric genotypes for all chromosomes
for(i in 1:7){
  chrom_idx = sites$CHROM == i
  chrom_geno = geno[,chrom_idx]
  write.table(t(chrom_geno),paste0("../Data/Genotype/B1K_red_chr_",i,"_num_mxn.txt"),quote = F)
}

# Extract the location of origin of the phenotyped individuals
phenotype_nonHSR = read.csv(phenotype_nonHSR_file)

genotype_location_map = phenotype_nonHSR %>%
  dplyr::select(Genotype,Location) %>%
  distinct(.keep_all = T)

K = read.csv("../Data/Genotype/B1K_GRM_red.csv", row.names = 1)

pca = prcomp(K)

# ---- Population structure analysis ----

# Create a data frame containing the loadings for the first 5 PCs 
df=data.frame(Genotype=rownames(K),
              Location=genotype_location_map$Location[match(rownames(K),genotype_location_map$Genotype)],
              PC1=pca$rotation[,1],
              PC2=pca$rotation[,2],
              PC3=pca$rotation[,3],
              PC4=pca$rotation[,4],
              PC5=pca$rotation[,5])

# Filter out the ambiguous accessions 
df_filter <- df[!is.na(df$Location),]

# Test homoscedasticity of locations in each PC 
leveneTest(PC1 ~ Location, data = df_filter) 
leveneTest(PC2 ~ Location, data = df_filter) # significant heteroscedasticity
leveneTest(PC3 ~ Location, data = df_filter) # significant heteroscedasticity
leveneTest(PC4 ~ Location, data = df_filter) # significant heteroscedasticity


# Test normality of standardized residuals in each PC 
test_stdres_normality(aov(PC1~Location,df_filter)) # significantly non-normal
test_stdres_normality(aov(PC2~Location,df_filter)) # significantly non-normal
test_stdres_normality(aov(PC3~Location,df_filter))
test_stdres_normality(aov(PC4~Location,df_filter)) # significantly non-normal

# Test separation of values due to location
kruskal.test(PC1~Location,data=df_filter) # significant separation
kruskal.test(PC2~Location,data=df_filter) # significant separation
kruskal.test(PC3~Location,data=df_filter) # significant separation
kruskal.test(PC4~Location,data=df_filter) # significant separation
kruskal.test(PC5~Location,data=df_filter) # significant separation


# ---- Population structure plot ---- 
df=data.frame(Genotype=rownames(K),Location=genotype_location_map$Location[match(rownames(K),genotype_location_map$Genotype)],PC1=pca$rotation[,1],PC2=pca$rotation[,2],PC3=pca$rotation[,3])

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

ggsave(paste0(figure_dir,"Kinship_Figure_new.png"),device = "png",bg = "white",width=12,height=5)

# ---- Linkage Disequilibrium ----
# The LD table was calculated MAF bin-wise 
# Computation of LD-decay distances was performed on a HPC cluster.
# Using: 
# LD_decay_HPC.R
# arg1: Input folder      - containing the TASSEL-generated LD tables with nomenclature: B1K_LD_chr[chr_nr]_full.txt
# arg2: Sample proportion - proportion [0-1] of marker pairs that should be sampled to compute function parameters
# arg3: Non-LD Distance   - distance where no LD is expected (~50Mb)
# arg4: Output folder     - folder where plots and analysis files are written.


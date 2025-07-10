## ---------------------------
##
## Script name: 5_GenomeMapping_GOenrichment
##
## Purpose of script:
##
## Author: M.Sc. Robin Lindner
##
## Date Created: 2025-07-10
##
## Copyright (c) Robin Lindner, 2025
## Email: robin.lindner@uni-potsdam.de
##
## ---------------------------
##
## Notes: Quick summary of genome mapping and GO enrichment pipeline.
##        
##        Variants were called on the outdated MorexV1 assembly.
##        => genome of MorexV1 and MorexV3 are aligned.
##        => SNP positions are adjusted to match MorexV3.
##   
##
## ---------------------------

## set working directory for Mac and PC

source("0_utils.R")

## ---------------------------

## ---- Genome alignment ----
# MorexV1 chromosome fasta files:
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/375/235/GCA_902375235.1_Morex_v1.0_update_x/GCA_902375235.1_Morex_v1.0_update_x_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/

# MorexV3 chromosome fasta files:
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/904/849/725/GCF_904849725.1_MorexV3_pseudomolecules_assembly/GCF_904849725.1_MorexV3_pseudomolecules_assembly_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/


# chromosomes were aligned using mummer:
# for(chrom in 1:7){
# nucmer --mum -c 100 -p ${chrom} MorexV3_Chromosomes/MorexV3_Seq_00${chrom}.fasta MorexV1_Chromosomes/MorexV1_Seq_00${chrom}.fasta
# }

# results in 7 coordinate files.

## ---- SNP mapping ----

coordinate_file_directory = ""

for(i in 1:7){
  chr_t = read.csv(file.path(coordinate_file_directory,"chrom",i,".coords.csv"),row.names=1)
  chr_SNP <- snps_MorexV1 %>%
    filter(Chromosome==i)
  res = map_SNP_positions(chr_t,chr_SNP)
  mappings = res %>%
    group_by(SNP) %>%
    filter(Seq_Identity == max(Seq_Identity)) %>%
    filter(!duplicated(SNP)) %>%
    mutate(Chromosome = i)
  if(i == 1){
    full_remap_df = mappings
  }else{
    full_remap_df = rbind(full_remap_df,mappings)
  }
}

write.csv(geno_remap_file,"../Data/Genotype/Assembly/SNP_remap.csv")



## ---- Gene ontology enrichment ----

GAF_file = "../Data/Genotype/GCF_904849725.1_MorexV3_pseudomolecules_assembly_gene_ontology.gaf"
GFF_file = "../Data/Genotype/Assembly/GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.gff"
Assembly_features_file = "../Data/Genotype/Assembly/GCF_904849725.1_MorexV3_pseudomolecules_assembly_feature_table.txt"

Annotation_outpath = "../Data/Genotype/Assembly/MorexV3_annotation.csv"


GO_df=readGAF(GAF_file)
MorexV3_anno = read.gff(GFF_file)

attr2df <- function(string){
  vec=unlist(strsplit(string,";"))
  named_list <- setNames(as.list(sub(".*=", "", vec)), sub("=.*", "", vec))
  df <- as.data.frame(named_list, stringsAsFactors = FALSE)
  return(df)
}
attr2df_v = Vectorize(attr2df)


anno_df=bind_rows(attr2df_v(MorexV3_anno$attributes))
full_df = cbind(MorexV3_anno,anno_df)
write.csv(full_df,Annotation_outpath,row.names = F)


assembly_features = read.table(Assembly_features_file,sep = "\t",header=F)
feature_classes=unlist(strsplit("feature	class	assembly	assembly_unit	seq_type	chromosome	genomic_accession	start	end	strand	product_accession	non-redundant_refseq	related_accession	name	symbol	GeneID	locus_tag	feature_interval_length	product_length	attributes",
                                split = "\\s+"))
names(assembly_features) = feature_classes

## 
ld_distances = c(12067,8746855,5617,4658,8110,451857,4944)

## ---- GO enrichment pSNP ---- 
## Define the relevant SNPs (top three)
snp_37_traits = "chr4H_632274504"

snp_29_traits = "chr5H_621066587"

snp_28_traits = "chr2H_704876309"

if(!dir.exists("../GO_enrichment_results")){
  dir.create("../GO_enrichment_results")
}
outdir="../GO_enrichment_results"

## 37 trait SNP
symbol_list_1=GeneIDs_for_cand_SNP(snp_37_traits,ld_distances,full_df)
res_table_1=FisherGOforAllDomains(symbol_list_1,geneList,geneID2GO,0.05)
prepareREVIGOIn(res_table_1,file.path(outdir,"GO_4_632Mb.txt"))

## 29 trait SNP
symbol_list_2=GeneIDs_for_cand_SNP(snp_29_traits,ld_distances,full_df)
res_table_2=FisherGOforAllDomains(symbol_list_2,geneList,geneID2GO,0.05)
prepareREVIGOIn(res_table_2,file.path(outdir,"GO_5_621Mb.txt"))


## 28 trait SNP
symbol_list_3=GeneIDs_for_cand_SNP(snp_28_traits,ld_distances,full_df)
res_table_3=FisherGOforAllDomains(symbol_list_3,geneList,geneID2GO,0.05)
prepareREVIGOIn(res_table_3,file.path(outdir,"GO_2_705Mb.txt"))

## ---- GO enrichment PH ----
ph_map = read.csv(ph_snp_ID_map)

## persistent
spec_snp_ID = c("ph-2-1","ph-2-3")
spec_snp = ph_map$SNP[ph_map$ID %in% spec_snp_ID]


spec_1=GeneIDs_for_cand_SNP(spec_snp[1],ld_distances,full_df,remap)
spec_2=GeneIDs_for_cand_SNP(spec_snp[2],ld_distances,full_df,remap)
res_table_spec1=FisherGOforAllDomains(spec_1,geneList,geneID2GO,0.05)
res_table_spec2=FisherGOforAllDomains(spec_2,geneList,geneID2GO,0.05)
prepareREVIGOIn(res_table_spec1,file.path(outdir,"GO_ph_2_1.txt"))
prepareREVIGOIn(res_table_spec2,file.path(outdir,"GO_ph_2_3.txt"))

## early
early_snp_ID = c("ph-1-1",
                 "ph-1-2",
                 "ph-1-3",
                 "ph-1-4",
                 "ph-1-5",
                 "ph-1-6",
                 "ph-1-7",
                 "ph-1-8",
                 "ph-1-9",
                 "ph-2-1",
                 "ph-4-1",
                 "ph-7-2")

early_snp = ph_map$SNP[ph_map$ID %in% early_snp_ID]


symbol_list_early=GeneIDs_for_cand_SNP(early_snp,ld_distances,full_df)
res_table_early=FisherGOforAllDomains(symbol_list_early,geneList,geneID2GO,0.05)
prepareREVIGOIn(res_table_early,file.path(outdir,"GO_ph_early.txt"))

## mid
mid_snp_ID = c("ph-1-10",
               "ph-2-1",
               "ph-2-3",
               "ph-3-2",
               "ph-7-3")

mid_snp = ph_map$SNP[ph_map$ID %in% mid_snp_ID]
symbol_list_mid=GeneIDs_for_cand_SNP(mid_snp,ld_distances,full_df,remap)
symbol_list_mid = na.omit(unique(symbol_list_mid))
res_table_mid=FisherGOforAllDomains(symbol_list_mid,geneList,geneID2GO,0.05)
prepareREVIGOIn(res_table_mid,file.path(outdir,"GO_ph_mid.txt"))


## late
late_snp_ID = c("ph-1-11",
                #"ph-2-1",
                "ph-2-2",
                #"ph-2-3",
                "ph-2-4",
                "ph-3-1",
                "ph-6-1",
                "ph-7-1",
                "ph-7-3")

late_snp = ph_map$SNP[ph_map$ID %in% late_snp_ID]
symbol_list_late=GeneIDs_for_cand_SNP(late_snp,ld_distances,full_df,remap)
symbol_list_late = na.omit(unique(symbol_list_late))
res_table_late=FisherGOforAllDomains(symbol_list_late,geneList,geneID2GO,0.05)
prepareREVIGOIn(res_table_late,file.path(outdir,"GO_ph_late.txt"))





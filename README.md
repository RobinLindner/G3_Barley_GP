# Barley phenome wide genomic association and HSR-supported GP

## Code
This directory contains all R scripts necessary to replicate the results from the publication.
### Main
0: file paths and functions used by other scripts. Is sourced at the beginning of all main scripts [1-5] <br>
1: Population structure and Linkage disequilibrium <br>
2: Genomic effec BLUPs and Heritability<br>
3: Phenome-wide time-series GWAS<br>
4: Genomic prediction of plant height, plant weight and NDVI score<br>
5: Genome mapping & Gene ontology enrichment for relevant markers (see main text)<br>

### Side
HPC_GWAS.R: allows user to performs GWAS in parallel on a high performance computing cluster.  <br>
HPC_ESA.R: allows user to extract significant associations at a set threshold from the GWAS results.<br>
MegaLMM_GP.R: allows user to run MegaLMM computations without marker fixed effects in parallel on HPC cluster.<br>
MegaLMM_MFE_GP.R: allows user to run MegaLMM computations with marker fixed effects in parallel on HPC cluster.<br>

## Data
This directory serves to structure the data as proposed in 0_utils.R. Genotyping-, phenotyping, as well as data generated in TASSEL will have to be added manually before scripts listed above will be able to be used.<br>
<br>
Alternatively, one can set file paths in 0_utils.R according to their own directory structure.<br>

## Supplements
trait_groups.csv:   linking traits to trait groups. <br>
ph_snp_map.csv:     linking PH SNPs to IDs used in the main text. <br>
GPgenotypes.txt:    subset of genotypes used for GP, necessary for parallel execution of MegaLMM. <br>
B1K_SNP_remap.csv:  mapping SNP positions of genotype from MorexV1 reference to MorexV3 (see main 5.). <br> 

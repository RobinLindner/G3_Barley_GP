# Barley phenome wide genomic association and HSR-supported GP

This repository represents an extensive documentation of computational analyses performed
for "Genetic architecture of temporal traits improves predictability of growth- and yield-related traits in wild barley."

## Folder structure:

###/Code
This directory contains all R scripts necessary to replicate the results from the publication.

Importantly, **"0_utils.R"** acts as a dictionary for file paths, in select cases 
these might have to be adjusted/consulted when incorporating data not stored on 
GitHub (i.e. Genotyping & Phenotyping data).

#### /1. Main results
This directory contains eight numbered scripts pertaining to different sections 
of the publication. These reproduce the main findings of the work, but might be
supported by supplementary scripts in other directories.


#### /2. Phenotype processing
This directory contains scripts that work in tandem with scripts "2_1_..." and 
"2_2_..." in /1.Main results. These scripts can be sourced locally or on a HPC 
cluster to decouple computation of BLUEs and BLUPs from the main results script.

#### /3. GWAS
This directory contains three scripts: <br>
HPC_GWAS.R: Perform a single GWAS with adjustable parameters. <br>
HPC_GWAS_exec.sh: Perform phenome-wide GWAS, in 1400+ scenarios in parallel on SLURM architecture. <br>
HPC_ESA.R: Extraction of Significant Associations; Parse the resulting GWAS files 
and extract significant marker trait associations, given a set significance threshold. <br>

####/4. Genomic prediction
This directory contains several scripts owing to the intricate CV scheme. 
All of the scripts utilize a CV mapping matrix, which splits the data into a 9:1
partition 100 times. This matrix was generated in script "4_GenomicPrediction.R", 
but upon replication of the results the provided matrix in "/Supplements/GP_CV_mat.csv"
should be used to ensure validity of the selected GP models. <br>

GWAS_Shuffle.R: Compute GWAS on a training set given the CV_mapping matrix, the run and the fold. <br>
GWAS_Shuffle_exec.R: Parallel execution of GWAS_Shuffle.R for all runs and folds. <br>
GWAS_Shuffle_eval.R: Evaluation of the GWAS results and identification of valid model scenarios. <br>

HPC_UVGP.R: Computes the univariate genomic predictions (GBLUP, HBLUP, G+HBLUP) for a given test set and scenario. <br>
UVGP_exec.sh: Parallel execution of HPC_UVGP for all valid scenarios and test sets. <br>
MegaLMM_GP.R: Predicts using MegaLMM without marker fixed effects for a given test set and scenario. <br>
MegaLMM_MFE_GP.R: Predicts using MegaLMM with marker fixed effects for a given test set and scenario. <br>
MVGP_exec.sh: Parallel execution of MegaLMM_GP.R and MegaLMM_MFE_GP.R for all valid scenarios and test sets. <br>

GP_postProcessing.R: Processing of GP results into tables and figures.<br>

####/5. LD decay
LD_decay.R: Estimation and plotting of LD decay for MAF bins. <br>

####/6. Morex V1->V3 genome alignment
genome_alignment.sh: alignment of full chromosome .fasta files using mummer (nucmer) <br>
Coord2CSV_cluster.R: transforms a given coordinate file resulting from genome_alignment.sh 
to a .csv format. <br>
Coord2CSV_cluster_exec.sh: parallel processing of all chromosomes via Coord2CSV_cluster.R <br>
SNPMapping.R: Use the fasta alignment to remap SNP to the aligned coordinates in MorexV3, selecting positions with maximal identity. <br>

####/7. Figures
The two scripts here were used to produce figures that were not plotted and saved 
in the main analysis scripts <br>

## Data
This directory serves to structure the data as proposed in 0_utils.R. Genotyping-, phenotyping, as well as data generated in TASSEL will have to be added manually before scripts listed above can be used.<br>
<br>
Alternatively, one can set file paths in 0_utils.R according to their own directory structure.<br>

## Supplements
trait_groups.csv:   linking traits to trait groups. <br>
ph_snp_map.csv:     linking PH SNPs to IDs used in the main text. <br>
GPgenotypes.txt:    subset of genotypes used for GP, necessary for parallel execution of MegaLMM. <br>
B1K_SNP_remap.csv:  mapping SNP positions of genotype from MorexV1 reference to MorexV3 (see main 5.). <br> 
GP_CV_mat_static.csv: the CV partitions that were used to compute the valid GP-models presented in the study.

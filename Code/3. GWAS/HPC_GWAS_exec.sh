#!/bin/bash
#SBATCH --job-name=BLUE_GWAS
#SBATCH --output=Debug/BLUE_GWAS.out
#SBATCH --error=Debug/BLUE_GWAS.err
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --mail-user=lindner5@uni-potsdam.de
#SBATCH --mail-type=FAIL

Rscript Code/HPC_GWAS.R Data/Generated/BLUEs_normalized.csv Data/Generated/GWAS_results Data/Genotype/B1K_red.vcf Data/Genotype/B1K_GRM_red.csv Data/Generated/GWAS_support
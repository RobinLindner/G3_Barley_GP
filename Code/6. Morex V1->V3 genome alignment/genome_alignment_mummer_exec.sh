#!/bin/bash
#SBATCH --job-name=genome_alignment           # Job name
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem-per-cpu=128G                    # Memory per node (specify how much memory you need per node)
#SBATCH --time=40:00:00             # Walltime (time limit)
#SBATCH --output=ga.out      # Standard output log file
#SBATCH --error=ga.err       # Standard error log file
#SBATCH --mail-user=lindner5@uni-potsdam.de
#SBATCH --mail-type=all

cd /work/lindner5/master/Thesis/Data/Gentoype/Assemblies
nucmer --maxmatch -c 50 -t 4 -p genome_alignment GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.fna GCA_902375235.1_Morex_v1.0_update_x_genomic.fna 

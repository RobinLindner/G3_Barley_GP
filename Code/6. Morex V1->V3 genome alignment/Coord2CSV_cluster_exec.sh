#!/bin/bash
#SBATCH --job-name=genome_alignment           # Job name
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --mem-per-cpu=8G                    # Memory per node (specify how much memory you need per node)
#SBATCH --time=40:00:00             # Walltime (time limit)
#SBATCH --output=Debug/c2c.out      # Standard output log file
#SBATCH --error=Debug/c2c.err       # Standard error log file
#SBATCH --mail-user=lindner5@uni-potsdam.de
#SBATCH --mail-type=all

cd /work/lindner5/master/Thesis/Data/Genotype/Assemblies

Rscript ../../../Code/Coord2CSV_cluster.R coord_files csv_files
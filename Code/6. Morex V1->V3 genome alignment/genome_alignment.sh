#!/bin/bash

# Define the strings and integers

# Nested loops
for chrom in {1..7}; do

  
  # Define a unique job name
  job_name="chromosome_${chrom}"
  
  # Optional: Print a log message for tracking
  echo "Submitting SLURM job: $job_name"

  # Create a SLURM job script dynamically
  cat <<EOT > ${job_name}.slurm
#!/bin/bash
#SBATCH --job-name=$job_name       # Job name
#SBATCH --output=Debug/${job_name}.out   # Standard output
#SBATCH --error=Debug/${job_name}.err    # Standard error
#SBATCH --time=08:00:00            # Time limit (HH:MM:SS)
#SBATCH --ntasks=1                 # Number of tasks
#SBATCH --cpus-per-task=2          # Number of CPU cores per task
#SBATCH --mem-per-cpu=80G                   # Memory per node (adjust as needed)

# Your actual code goes here
cd /work/lindner5/master/Thesis/Data/Genotype/Assemblies
nucmer --mum -c 100 -p ${chrom} MorexV3_Chromosomes/MorexV3_Seq_00${chrom}.fasta MorexV1_Chromosomes/MorexV1_Seq_00${chrom}.fasta
wait
# Replace this with the actual command(s) to run
EOT

  # Submit the SLURM job
  sbatch ${job_name}.slurm

  # Clean up generated job script 
  rm ${job_name}.slurm

done
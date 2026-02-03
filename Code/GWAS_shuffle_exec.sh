#!/bin/bash

# Define the strings and integers
traits=("RGB1_Plant_Avg_HEIGHT_MM" "VNIR_Plant_NDVI.avg" "SC_Plant_Weight")  # Replace with your strings
dats=(14 21 28 35 42)                        # Replace with your integers

# Nested loops
for trait in "${traits[@]}"; do
for dat in "${dats[@]}"; do
for i in {1..100}; do

output_file="Data/Generated/GenomicPrediction/GWAS_Shuffle/${trait}_dat${dat}_run${i}_sigSNP.csv"
if [ -f "$output_file1" ] ; then
echo "Files $output_file1 already exist. Skipping iteration."
continue
fi

# Define a unique job name
job_name="run_${trait}_${dat}_${i}"

# Optional: Print a log message for tracking
echo "Submitting SLURM job: $job_name"

mamba activate G3_pub

# Create a SLURM job script dynamically
cat <<EOT > ${job_name}.slurm
#!/bin/bash
#SBATCH --job-name=$job_name       # Job name
#SBATCH --output=Debug/${job_name}.out   # Standard output
#SBATCH --error=Debug/${job_name}.err    # Standard error
#SBATCH --time=03:00:00            # Time limit (HH:MM:SS)
#SBATCH --ntasks=2                 # Number of tasks
#SBATCH --cpus-per-task=2          # Number of CPU cores per task
#SBATCH --mem-per-cpu=6G                   # Memory per cpu (adjust as needed)
#SBATCH --mail-user=lindner5@uni-potsdam.de
#SBATCH --mail-type=FAIL

# Your actual code goes here
cd /work/lindner5/master/G3_Barley_GP

srun --ntasks=1 Rscript Code/GWAS_Shuffle Data/Supplements/GP_CV_mat2.csv Data/Generated/GenomicPrediction/GWAS_Shuffle ${trait} ${dat} ${i}  &
  wait
  EOT

# Submit the SLURM job
sbatch ${job_name}.slurm

# Clean up generated job script 
rm ${job_name}.slurm

# rm Data/Generated/GenomicPrediction/Results/Runs/${trait}_dat${dat}_run${i}_MegaLMM_fold*
# rm Data/Generated/GenomicPrediction/Results_MFE/Runs/${trait}_dat${dat}_run${i}_MegaLMM_fold*

done
done
sleep 1h
done
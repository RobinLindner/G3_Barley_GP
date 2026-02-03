#!/bin/bash
#SBATCH --job-name=GWAS_Shuffle
#SBATCH --output=Debug/run_%A_%a.out
#SBATCH --error=Debug/run_%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --mail-user=lindner5@uni-potsdam.de
#SBATCH --mail-type=FAIL
#SBATCH --array=0-1499%20          # 1500 jobs total, 20 at a time

# 1. Define your arrays
traits=("RGB1_Plant_Avg_HEIGHT_MM" "VNIR_Plant_NDVI.avg" "SC_Plant_Weight")
dats=(14 21 28 35 42)
num_runs=100

# 2. Map the Task ID to your variables using math
# This replaces your nested loops
trait_idx=$(( SLURM_ARRAY_TASK_ID / (5 * num_runs) ))
remainder=$(( SLURM_ARRAY_TASK_ID % (5 * num_runs) ))
dat_idx=$(( remainder / num_runs ))
run_num=$(( (remainder % num_runs) + 1 ))

trait=${traits[$trait_idx]}
dat=${dats[$dat_idx]}

# 3. Skip if file exists
output_file="Data/Generated/GenomicPrediction/GWAS_Shuffle/${trait}_dat${dat}_run${run_num}_sigSNP.csv"
if [ -f "$output_file" ]; then
echo "File $output_file exists. Skipping."
exit 0
fi

# 4. Execution
# Note: Source your conda/mamba profile if 'mamba activate' isn't in your .bashrc
source $(conda info --base)/etc/profile.d/conda.sh
conda activate G3_pub

cd /work/lindner5/master/G3_Barley_GP

echo "Processing Trait: $trait | Dat: $dat | Run: $run_num"
Rscript Code/GWAS_Shuffle Data/Supplements/GP_CV_mat2.csv Data/Generated/GenomicPrediction/GWAS_Shuffle ${trait} ${dat} ${run_num}

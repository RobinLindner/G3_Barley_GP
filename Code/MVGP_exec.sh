#!/bin/bash
#SBATCH --job-name=MultiVariateGP
#SBATCH --output=Debug/run_%A_%a.out
#SBATCH --error=Debug/run_%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --mail-user=lindner5@uni-potsdam.de
#SBATCH --mail-type=FAIL


# 1. Define your arrays
traits=("RGB1_Plant_Avg_HEIGHT_MM" "VNIR_Plant_NDVI.avg" "SC_Plant_Weight")
dats=(14 21 28 35 42)


# 2. Map the Task ID to your variables using math
# This replaces your nested loops

ACTUAL_ID=$(( SLURM_ARRAY_TASK_ID + ${OFFSET:-0} ))

trait_idx=$(( ACTUAL_ID / (5) ))
remainder=$(( ACTUAL_ID % (5) ))
dat_idx=$(( remainder ))

scenarioID=$(( ACTUAL_ID + 1 ))
trait=${traits[$trait_idx]}
dat=${dats[$dat_idx]}

# 3. Skip if file exists

output_path="Data/Generated/GenomicPrediction/MVGP/"
output_file="Data/Generated/GenomicPrediction/MVGP/Accuracy/${trait}_${dat}_accuracy.csv"
if [ -f "$output_file" ]; then
echo "File $output_file exists. Skipping."
exit 0
fi

# 4. Execution
# Note: Source your conda/mamba profile if 'mamba activate' isn't in your .bashrc
#source $(conda info --base)/etc/profile.d/conda.sh
#conda activate G3_pub

cd /work/lindner5/master/G3_Barley_GP

echo "Processing Trait: $trait | Dat: $dat"
Rscript Code/MegaLMM_GP.R ${output_path} ${trait} ${dat} $scenarioID
Rscript Code/MegaLMM_MFE_GP.R ${output_path} ${trait} ${dat} $scenarioID
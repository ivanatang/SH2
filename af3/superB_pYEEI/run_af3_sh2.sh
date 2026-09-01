#!/bin/bash
#SBATCH --job-name=af3_sh2_pYEEI
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=02:30:00
#SBATCH --partition=al40
#SBATCH --qos=gpu-normal
#SBATCH --gres=gpu:l40:1
#SBATCH --account=ucb351_asc4
#SBATCH --output=output_af3_sh2_%j.out
#SBATCH --error=error_af3_sh2_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ivana.tang@colorado.edu

module purge
module load alphafold/3.0.0

export INPUT_FILE=$(pwd)/superB_pYEEI/superbinder_sh2_pYEEI.json
export OUTPUT_DIR=$(pwd)/superB_pYEEI/output
export AF3_MODEL_PARAMETERS_DIR=/projects/$USER/alphafold3_weights/   # path where you stored the downloaded AF3 weights

mkdir -p "$OUTPUT_DIR"

run_alphafold \
  --json_path=$INPUT_FILE \
  --output_dir=$OUTPUT_DIR \
  --model_dir=$AF3_MODEL_PARAMETERS_DIR

#!/bin/bash
#SBATCH --partition=aa100
#SBATCH --gres=gpu:a100-40gb:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --qos=gpu-normal
#SBATCH --time=02:00:00
#SBATCH --job-name=boltz_cSrc_validation
#SBATCH --output=logs/boltz_cSrc_%j.out
#SBATCH --error=logs/boltz_cSrc_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ivta1597@colorado.edu

module purge
module load anaconda
conda activate boltz2

boltz predict /projects/$USER/SH2/boltz2/WT_SH2_pYEEI/cSrc_SH2_EPQpYEEIPIYL.yaml \
    --out_dir /projects/$USER/SH2/boltz2/WT_SH2_pYEEI/boltz_cSrc_validation \
    --use_msa_server \
    --use_potentials \
    --diffusion_samples 25 \
    --recycling_steps 10 \
    --seed 1 \
    --no_kernels

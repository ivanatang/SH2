#!/bin/bash
#SBATCH --partition=aa100
#SBATCH --gres=gpu:a100_3g.20gb:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --qos=gpu-testing
#SBATCH --time=00:15:00
#SBATCH --job-name=boltz_test
#SBATCH --output=boltz_test.%j.out

module purge
module load anaconda
conda activate boltz2

boltz predict /projects/$USER/SH2/boltz2/cSrc_SH2_EPQpYEEIPIYL.yaml \
    --out_dir /scratch/alpine/$USER/boltz_test \
    --use_msa_server \
    --use_potentials \
    --diffusion_samples 1 \
    --recycling_steps 3 \
    --seed 1 \
    --no_kernels

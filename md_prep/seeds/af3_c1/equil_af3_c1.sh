#!/bin/bash

#SBATCH --job-name=eq_af3_c1
#SBATCH --output=output_eq_af3_c1_%j.out
#SBATCH --error=error_eq_af3_c1_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

export TMPDIR=$SLURM_SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load gcc
module load openmpi
module load anaconda
module load gromacs

conda activate SH2

DIR=/projects/ivta1597/SH2/md_prep
MDP=$DIR/mdp
TOP=$DIR/seeds/af3_c1/gromacs/af3_c1_dodecahedron.top
SEQ_DIR=/scratch/alpine/ivta1597/SH2/seeds/af3_c1

# NVT
cd "$SEQ_DIR"
mkdir -p NVT
cd NVT
gmx grompp -f $MDP/nvt.mdp -c $SEQ_DIR/EM/em.gro -r $SEQ_DIR/EM/em.gro -p $TOP -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT
cd "$SEQ_DIR"
mkdir -p NPT
cd NPT
gmx grompp -f $MDP/npt.mdp -c $SEQ_DIR/NVT/nvt.gro -t $SEQ_DIR/NVT/nvt.cpt -r $SEQ_DIR/NVT/nvt.gro -p $TOP -o npt.tpr
gmx mdrun -deffnm npt

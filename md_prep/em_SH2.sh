#!/bin/bash

#SBATCH --job-name=em_SH2
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
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
SEQ_DIR=/scratch/alpine/ivta1597/SH2

mkdir -p "$SEQ_DIR/EM"
cd "$SEQ_DIR/EM"
gmx grompp -f $MDP/em.mdp \
  -c $DIR/gromacs/sh2_superbinder_pTyr_dodecahedron.gro \
  -p $DIR/gromacs/sh2_superbinder_pTyr_dodecahedron.top \
  -o em.tpr
gmx mdrun -deffnm em

#!/bin/bash

#SBATCH --job-name=SH2_prod
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
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
TOP=$DIR/gromacs/sh2_superbinder_pTyr_dodecahedron.top
SEQ_DIR=/scratch/alpine/ivta1597/SH2

cd "$SEQ_DIR"
mkdir -p prod_md_1x${SLURM_CPUS_PER_TASK}
cd prod_md_1x${SLURM_CPUS_PER_TASK}
gmx grompp -f $MDP/prod_md.mdp -c $SEQ_DIR/NPT/npt.gro -t $SEQ_DIR/NPT/npt.cpt -p $TOP -o prod_md_500ns.tpr

# Single-node thread-MPI (no mpirun/gmx_mpi), matching the already-working
# em_SH2.sh/equil_SH2.sh pattern -- NOT the real-MPI (gmx_mpi + mpirun -np 64)
# approach borrowed from prod_md_PYR1_LCA.sh's 64-rank / hand-tuned DD grid.
# That approach caused a real, observed failure: the first attempt at this run
# only reached step 110720 (221 ps of the requested 50000 ps) after most of a
# 12-hour allocation -- ~0.46 ns/day, with the .log showing severe DD load
# imbalance ("force 143.3%") and a bad PME/force ratio. 64 real MPI ranks on a
# single node is very likely over-decomposed for a ~21k-atom system; this
# job was always --nodes=1 anyway, so real inter-process MPI was pure overhead
# for zero benefit. No explicit -ntomp/-ntmpi here, same as em/equil -- relies
# on OMP_NUM_THREADS above and thread-MPI's own auto-tuning, which is the
# GROMACS-recommended default for single-node runs (unlike a hand-guessed DD
# grid for real MPI, this auto-tuning is specifically designed for this case).
gmx mdrun -deffnm prod_md_500ns

#!/bin/bash

#SBATCH --job-name=SH2_prod
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --constraint=ib
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

export TMPDIR=$SLURM_SCRATCH
export SLURM_EXPORT_ENV=ALL

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
mkdir -p prod_md_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}
cd prod_md_${SLURM_NTASKS}x${SLURM_CPUS_PER_TASK}
gmx_mpi grompp -f $MDP/prod_md.mdp -c $SEQ_DIR/NPT/npt.gro -t $SEQ_DIR/NPT/npt.cpt -p $TOP -o prod_md_500ns.tpr

# NOTE: unlike prod_md_PYR1_LCA.sh, this deliberately does NOT pass -dd/-npme.
# Their script hand-tunes those (DD grid 6x4x2 + 16 PME ranks) for the PYR1+LCA
# system's specific size/box -- I have no way to verify that same grid is valid
# for this (smaller, differently-shaped dodecahedron) system without a real
# MPI-enabled GROMACS build to test against (this Mac's local `gmx` is
# thread-MPI only, no true multi-rank MPI). An invalid manual DD grid makes
# mdrun fail immediately, which you do not want discovered 20 hours into a
# 24-hour allocation. Letting mdrun auto-decompose is slower to start (it logs
# its chosen grid near the top of the .log) but guaranteed valid. Once you've
# seen what it picks and how it scales, this can be replaced with an explicit,
# actually-verified -dd/-npme for future runs.
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK

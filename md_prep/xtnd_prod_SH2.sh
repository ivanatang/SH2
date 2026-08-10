#!/bin/bash

#SBATCH --job-name=SH2_prod_xtnd
#SBATCH --output=output_xtnd_%j.out
#SBATCH --error=error_xtnd_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#
# Resumes an interrupted prod_md_SH2.sh run from its checkpoint, continuing
# toward the SAME original target (prod_md.mdp's nsteps=25000000, i.e. 50 ns --
# this is NOT extending the run past 50 ns, just finishing the run that ran
# out of wall time before reaching it). Based on biosensors' xtnd_prod_PYR1_LCA.sh
# (-s / -cpi / -append pattern), but using this project's actual benchmark-
# determined settings (48-rank real MPI, ntomp=1, auto DD/PME -- see
# prod_md_SH2.sh) rather than their hand-tuned -npme/-dd, which was never
# benchmarked for this system and is exactly what caused the earlier 0.46 ns/day
# failure at a different rank count.
#
# Usage: sbatch xtnd_prod_SH2.sh
# (resubmit again, unmodified, if it runs out of wall time again before 50 ns)

export TMPDIR=$SLURM_SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load gcc
module load openmpi
module load anaconda
module load gromacs

conda activate SH2

SEQ_DIR=/scratch/alpine/ivta1597/SH2

cd "$SEQ_DIR/prod_md_48x1"

if [ ! -f prod_md_500ns.cpt ]; then
    echo "ERROR: no checkpoint found at $SEQ_DIR/prod_md_48x1/prod_md_500ns.cpt -- nothing to resume from." >&2
    exit 1
fi

mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns \
  -s prod_md_500ns.tpr -cpi prod_md_500ns.cpt -append -ntomp $SLURM_CPUS_PER_TASK

#!/bin/bash

#SBATCH --job-name=boltz_c1_prod_xtnd
#SBATCH --output=output_resume_boltz_c1_%j.out
#SBATCH --error=error_resume_boltz_c1_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#
# Resumes an interrupted prod_md_boltz_c1.sh run from its checkpoint, continuing
# toward the SAME original target (prod_md_200ns.mdp's nsteps=100000000, i.e.
# 200 ns -- NOT extending past 200 ns, just finishing the run that ran out of
# wall time). Mirrors xtnd_prod_SH2.sh exactly (same -s/-cpi/-append pattern,
# same --exclusive -- added there after a resumed run measured ~0.35 ns/day
# from node-sharing DD load imbalance; kept here from the start).
#
# Usage: sbatch resume_prod_boltz_c1.sh
# (resubmit again, unmodified, if it runs out of wall time again before 200 ns)

export TMPDIR=$SLURM_SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load gcc
module load openmpi
module load anaconda
module load gromacs

conda activate SH2

SEQ_DIR=/scratch/alpine/ivta1597/SH2/seeds/boltz_c1

cd "$SEQ_DIR/prod_md_48x1"

if [ ! -f prod_md_200ns.cpt ]; then
    echo "ERROR: no checkpoint found at $SEQ_DIR/prod_md_48x1/prod_md_200ns.cpt -- nothing to resume from." >&2
    exit 1
fi

mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_200ns \
  -s prod_md_200ns.tpr -cpi prod_md_200ns.cpt -append -ntomp $SLURM_CPUS_PER_TASK

#!/bin/bash

#SBATCH --job-name=af3_c2_prod_xtnd
#SBATCH --output=output_resume_af3_c2_%j.out
#SBATCH --error=error_resume_af3_c2_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=22:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#
# Resumes an interrupted prod_md_af3_c2.sh run from its checkpoint, continuing
# toward the SAME original target (prod_md_200ns.mdp's nsteps=100000000, i.e.
# 200 ns -- NOT extending past 200 ns, just finishing the run that ran out of
# wall time). Mirrors xtnd_prod_SH2.sh's -s/-cpi/-append pattern, but WITHOUT
# --exclusive (matching prod_md_af3_c2.sh) -- this partition/qos rejects it
# at sbatch time ("Error 22: The oversubscribe Slurm directive cannot be
# used with the provided partition"), unlike the original single-system
# scripts where it worked. xtnd_prod_SH2.sh's own comment documents a
# ~0.35 ns/day node-sharing slowdown without --exclusive, so check actual
# throughput with check_prod_performance.sh rather than assuming it's fine.
#
# Usage: sbatch resume_prod_af3_c2.sh
# (resubmit again, unmodified, if it runs out of wall time again before 200 ns)

export TMPDIR=$SLURM_SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# gcc/openmpi/anaconda modules were removed from Alpine's module stack
# (confirmed via `module spider` on 2026-09-04 -- only "gromacs" remains).
# Point directly at the real installs instead of relying on modules for these three
# (gromacs.lua's own internal load('gcc')/load('openmpi') calls also now fail, but
# harmlessly -- its prepend_path calls for gmx_mpi still run regardless).
export PATH=/curc/sw/install/openmpi/5.0.6/gcc/14.2.0/bin:$PATH
export LD_LIBRARY_PATH=/curc/sw/install/gcc/14.2.0/lib64:/curc/sw/install/openmpi/5.0.6/gcc/14.2.0/lib:$LD_LIBRARY_PATH
source /curc/sw/install/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate SH2

module purge
module load gromacs

SEQ_DIR=/scratch/alpine/ivta1597/SH2/seeds/af3_c2

cd "$SEQ_DIR/prod_md_48x1"

if [ ! -f prod_md_200ns.cpt ]; then
    echo "ERROR: no checkpoint found at $SEQ_DIR/prod_md_48x1/prod_md_200ns.cpt -- nothing to resume from." >&2
    exit 1
fi

mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_200ns \
  -s prod_md_200ns.tpr -cpi prod_md_200ns.cpt -append -ntomp $SLURM_CPUS_PER_TASK

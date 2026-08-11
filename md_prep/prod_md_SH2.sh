#!/bin/bash

#SBATCH --job-name=SH2_prod
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err
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
mkdir -p prod_md_48x1
cd prod_md_48x1
gmx grompp -f $MDP/prod_md.mdp -c $SEQ_DIR/NPT/npt.gro -t $SEQ_DIR/NPT/npt.cpt -p $TOP -o prod_md_500ns.tpr

# Real MPI, 48 ranks, ntomp=1, DD/PME auto-tuned (no manual -dd/-npme) --
# determined by the benchmark sweep in bchmk_prod_md_SH2.sh/submit_benchmarks_SH2.sh
# (see parse_bchmk_performance_SH2.py output): 48-rank real MPI measured at
# 260.36 +/- 2.54 ns/day (mean of 2 reps), clearly best across the full sweep
# (4-64 cores, both real MPI and single-node thread-MPI). Notably, 64 ranks
# measured WORSE (146.38 ns/day) than 48 -- consistent with the earlier
# catastrophic 64-rank failure (~0.46 ns/day) on the first production attempt,
# which is why this was benchmarked rather than assumed. Single-node thread-MPI
# (previously used here as a safe fallback after that failure) topped out at
# 170.25 ns/day at 64 cores -- worse than 48-rank real MPI at every point
# measured, so this is not simply reverting to what failed before; it's the
# actual measured optimum, not the same untested guess that failed originally.
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_500ns -ntomp $SLURM_CPUS_PER_TASK

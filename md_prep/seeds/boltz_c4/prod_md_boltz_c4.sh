#!/bin/bash

#SBATCH --job-name=boltz_c4_prod
#SBATCH --output=output_prod_boltz_c4_%j.out
#SBATCH --error=error_prod_boltz_c4_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=22:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

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

DIR=/projects/ivta1597/SH2/md_prep
MDP=$DIR/mdp
TOP=$DIR/seeds/boltz_c4/gromacs/boltz_c4_dodecahedron.top
SEQ_DIR=/scratch/alpine/ivta1597/SH2/seeds/boltz_c4

cd "$SEQ_DIR"
mkdir -p prod_md_48x1
cd prod_md_48x1
gmx grompp -f $MDP/prod_md_200ns.mdp -c $SEQ_DIR/NPT/npt.gro -t $SEQ_DIR/NPT/npt.cpt -p $TOP -o prod_md_200ns.tpr

# Same benchmark-determined settings as prod_md_SH2.sh: 48-rank real MPI,
# ntomp=1, auto DD/PME (measured 260.36 +/- 2.54 ns/day for this system size
# on this partition -- see prod_md_SH2.sh's comment for the full benchmark
# rationale). NOTE: unlike prod_md_SH2.sh/xtnd_prod_SH2.sh, --exclusive is
# NOT used here -- sbatch rejects it on this partition/qos ("Error 22: The
# oversubscribe Slurm directive cannot be used with the provided partition").
# xtnd_prod_SH2.sh's own comment says --exclusive was added after a run
# without it measured ~0.35 ns/day from node-sharing DD load imbalance, so
# throughput here isn't guaranteed to match the 260.36 ns/day benchmark --
# check it with check_prod_performance.sh once running, don't just assume.
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_200ns -ntomp $SLURM_CPUS_PER_TASK

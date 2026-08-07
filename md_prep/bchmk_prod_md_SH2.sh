#!/bin/bash

#SBATCH --job-name=bchmk_SH2
#SBATCH --output=output_benchmark_%j.out
#SBATCH --error=error_benchmark_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=FAIL
#
# --ntasks / --cpus-per-task are intentionally NOT set here -- submit_benchmarks_SH2.sh
# overrides them per sweep point at submission time (`sbatch --ntasks=N --cpus-per-task=M ...`).
#
# Usage (normally called by submit_benchmarks_SH2.sh, not run directly):
#   sbatch --ntasks=N --cpus-per-task=M bchmk_prod_md_SH2.sh <mode> <rep>
#     mode = "mpi"    -> real MPI (gmx_mpi + mpirun -np N), ntomp=1, auto DD/PME
#            "thread" -> single-node thread-MPI (plain gmx, no mpirun), ntomp=M
#     rep  = replicate number, just used in the output directory name

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

MODE=$1   # mpi | thread
REP=$2

if [ "$MODE" != "mpi" ] && [ "$MODE" != "thread" ]; then
    echo "ERROR: mode must be 'mpi' or 'thread', got '$MODE'" >&2
    exit 1
fi

if [ "$MODE" == "mpi" ]; then
    CORES=$SLURM_NTASKS
else
    CORES=$SLURM_CPUS_PER_TASK
fi

OUTDIR=$SEQ_DIR/benchmark/${MODE}_${CORES}cores_rep${REP}
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Same NPT starting point every time so all sweep points are directly comparable.
gmx grompp -f $MDP/prod_md_benchmark.mdp -c $SEQ_DIR/NPT/npt.gro -t $SEQ_DIR/NPT/npt.cpt -p $TOP -o prod_bench.tpr

if [ "$MODE" == "mpi" ]; then
    # Real MPI, ntomp=1 per rank, DD/PME left to auto-tune (no -dd/-npme) --
    # matches the biosensors benchmarking methodology, but at core counts
    # appropriate to our system's actual size (see submit_benchmarks_SH2.sh).
    mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_bench -ntomp 1
else
    # Single-node thread-MPI, no mpirun/gmx_mpi -- same pattern as the (fixed)
    # prod_md_SH2.sh. No explicit -ntmpi/-ntomp: let thread-MPI auto-tune the
    # rank/thread split for the requested core count.
    gmx mdrun -deffnm prod_bench
fi

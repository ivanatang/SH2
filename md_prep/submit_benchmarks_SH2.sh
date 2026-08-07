#!/bin/bash
# submit_benchmarks_SH2.sh
# ─────────────────────────────────────────────────────────────────────────────
# Submits the full core/thread benchmark sweep for bchmk_prod_md_SH2.sh.
#
# Tests two axes at each core count in CORE_COUNTS:
#   mpi    -- real MPI (gmx_mpi + mpirun -np N), ntomp=1, DD/PME auto-tuned
#   thread -- single-node thread-MPI (plain gmx), no mpirun
#
# Core counts are weighted toward the low end on purpose: this system is
# ~21k atoms, roughly half the size of the biosensors PYR1+LCA systems whose
# own benchmarking found 64-rank real MPI worked well -- for our system, the
# first real production attempt showed 64-rank real MPI performing
# catastrophically (~0.46 ns/day, severe DD load imbalance), almost certainly
# because domains that fine-grained are too small for a system this size. The
# actual optimum is very plausibly well below 64 cores; this sweep is designed
# to find it rather than assume it.
#
# Usage:
#   bash submit_benchmarks_SH2.sh [n_reps]     # n_reps defaults to 2
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

N_REPS=${1:-2}
CORE_COUNTS=(4 8 12 16 24 32 48 64)
SCRIPT=bchmk_prod_md_SH2.sh

if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: $SCRIPT not found in $(pwd) -- run this from /projects/ivta1597/SH2/md_prep" >&2
    exit 1
fi

echo "Submitting benchmark sweep: ${#CORE_COUNTS[@]} core counts x 2 modes x ${N_REPS} reps = $(( ${#CORE_COUNTS[@]} * 2 * N_REPS )) jobs"

for cores in "${CORE_COUNTS[@]}"; do
    for rep in $(seq 1 "$N_REPS"); do
        # Real MPI: N ranks, ntomp=1 each
        jid_mpi=$(sbatch --parsable --ntasks="$cores" --cpus-per-task=1 \
            --job-name="bench_mpi_${cores}" \
            "$SCRIPT" mpi "$rep")
        echo "  mpi    ${cores} cores rep${rep}: job $jid_mpi"

        # Thread-MPI: 1 task, N OpenMP/thread-MPI threads
        jid_thread=$(sbatch --parsable --ntasks=1 --cpus-per-task="$cores" \
            --job-name="bench_thread_${cores}" \
            "$SCRIPT" thread "$rep")
        echo "  thread ${cores} cores rep${rep}: job $jid_thread"
    done
done

echo ""
echo "All jobs submitted. Check status with:"
echo "  squeue -u \$USER --format=\"%.10i %.20j %.10T %.20R\""
echo ""
echo "Once they finish, parse results with:"
echo "  python3 parse_bchmk_performance_SH2.py"

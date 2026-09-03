#!/bin/bash
#
# Check progress and throughput of ongoing/finished production runs for all
# 8 seeds. Run this ON THE CLUSTER (uses /scratch/alpine paths, gmx, and
# sacct). Companion to check_em_eq_status.sh -- use that one first for
# EM/equilibration, this one once production is running.
#
# Why this matters here specifically: prod_md_<seed>.sh and
# resume_prod_<seed>.sh do NOT use --exclusive (this partition/qos rejects
# it -- "Error 22: The oversubscribe Slurm directive cannot be used with
# the provided partition"), unlike the original single-system
# prod_md_SH2.sh/xtnd_prod_SH2.sh. xtnd_prod_SH2.sh's own comment documents
# a run without --exclusive measuring ~0.35 ns/day (vs. the 260.36 ns/day
# benchmark) from node-sharing DD load imbalance -- so actual throughput
# here needs to be checked, not assumed to match the benchmark.
#
# For each seed, reports:
#   - SLURM job state via sacct (most recent job named <seed>_prod or
#     <seed>_prod_xtnd, whichever is more recent)
#   - current step / progress toward the 200 ns (100000000 step) target,
#     from the most recent "Writing checkpoint, step N" line in the log
#     (checkpointed periodically during the run, so this works even for
#     a still-running job -- no need to wait for completion)
#   - realized throughput: current_ns / elapsed_wall_time, compared against
#     the 260.36 ns/day benchmark -- flags anything running much slower
#     (likely node-sharing contention, exactly the risk noted above)
#   - estimated time remaining at the realized rate
#   - Fatal error check
#
# Usage: bash check_prod_performance.sh [seed ...]
#   (no args = check all 8; or pass specific seed names, e.g. af3_c1 boltz_c2)

module purge
module load gcc
module load openmpi
module load anaconda
module load gromacs
conda activate SH2

SCRATCH_ROOT=/scratch/alpine/ivta1597/SH2/seeds
ALL_SEEDS="af3_c1 af3_c2 af3_c3 af3_c4 boltz_c1 boltz_c2 boltz_c3 boltz_c4"
SEEDS="${*:-$ALL_SEEDS}"

TARGET_STEPS=100000000   # prod_md_200ns.mdp nsteps
DT_NS=0.000002           # 0.002 ps timestep, in ns
BENCHMARK_NS_DAY=260.36  # prod_md_SH2.sh benchmark, --exclusive, for comparison only

sacct_job () {
  # $1 = job name pattern (matches "<seed>_prod" or "<seed>_prod_xtnd") --
  # prints "State Elapsed Start" for the most recent matching job, tab-separated
  sacct -u "$USER" --name="$1" --format=State,Elapsed,Start --noheader -X 2>/dev/null | tail -1
}

elapsed_to_hours () {
  # SLURM Elapsed is [D-]HH:MM:SS -- convert to decimal hours
  local e="$1" days=0 rest="$1"
  if [[ "$e" == *-* ]]; then
    days="${e%%-*}"
    rest="${e#*-}"
  fi
  local h m s
  IFS=: read -r h m s <<< "$rest"
  echo "scale=4; $days*24 + ${h#0} + ${m#0}/60 + ${s#0}/3600" | bc 2>/dev/null
}

printf "%-10s | %-10s %-10s %-12s %-8s | %-12s %-9s %-10s\n" \
  "seed" "job" "state" "elapsed(h)" "progress" "current_ns" "ns/day" "ETA(h)"
printf '%.0s-' {1..100}; echo

for seed in $SEEDS; do
  DIR="$SCRATCH_ROOT/$seed/prod_md_48x1"
  log="$DIR/prod_md_200ns.log"

  # prefer the resume job if one exists (most recent attempt), else the initial prod job
  resume_info=$(sacct_job "${seed}_prod_xtnd")
  prod_info=$(sacct_job "${seed}_prod")
  if [ -n "$resume_info" ]; then
    job_label="resume"; info="$resume_info"
  else
    job_label="prod"; info="$prod_info"
  fi

  if [ -z "$info" ]; then
    printf "%-10s | %-10s %-10s %-12s %-8s | %-12s %-9s %-10s\n" "$seed" "$job_label" "no job found" "-" "-" "-" "-" "-"
    continue
  fi

  state=$(echo "$info" | awk '{print $1}')
  elapsed=$(echo "$info" | awk '{print $2}')
  elapsed_h=$(elapsed_to_hours "$elapsed")

  if [ ! -f "$log" ]; then
    printf "%-10s | %-10s %-10s %-12s %-8s | %-12s %-9s %-10s\n" "$seed" "$job_label" "$state" "$elapsed_h" "no log yet" "-" "-" "-"
    continue
  fi

  if grep -qi "Fatal error" "$log"; then
    printf "%-10s | %-10s %-10s %-12s %-8s | %-12s %-9s %-10s\n" "$seed" "$job_label" "$state" "$elapsed_h" "FATAL_ERROR" "-" "-" "-"
    continue
  fi

  last_step=$(grep -oE "Writing checkpoint, step [0-9]+" "$log" | tail -1 | awk '{print $NF}')
  if [ -z "$last_step" ]; then
    # no checkpoint written yet (run just started, or still in grompp) --
    # fall back to the final summary if the run already finished cleanly
    last_step=$(grep -A1 "^ *Step *Time" "$log" | tail -1 | awk '{print $1}')
  fi

  if [ -z "$last_step" ]; then
    printf "%-10s | %-10s %-10s %-12s %-8s | %-12s %-9s %-10s\n" "$seed" "$job_label" "$state" "$elapsed_h" "no progress yet" "-" "-" "-"
    continue
  fi

  current_ns=$(echo "scale=3; $last_step * $DT_NS" | bc)
  pct=$(echo "scale=1; $last_step * 100 / $TARGET_STEPS" | bc)

  ns_day="N/A"; eta_h="N/A"
  if [ -n "$elapsed_h" ] && [ "$(echo "$elapsed_h > 0" | bc)" = "1" ]; then
    ns_day=$(echo "scale=2; $current_ns / ($elapsed_h / 24)" | bc)
    remaining_ns=$(echo "scale=3; ($TARGET_STEPS - $last_step) * $DT_NS" | bc)
    if [ "$(echo "$ns_day > 0" | bc)" = "1" ]; then
      eta_h=$(echo "scale=1; $remaining_ns / $ns_day * 24" | bc)
    fi
    slow_flag=""
    if [ "$(echo "$ns_day < $BENCHMARK_NS_DAY / 2" | bc)" = "1" ]; then
      slow_flag=" <-- SLOW (< 50% of ${BENCHMARK_NS_DAY} ns/day benchmark; check for node-sharing contention)"
    fi
  fi

  printf "%-10s | %-10s %-10s %-12s %-8s | %-12s %-9s %-10s%s\n" \
    "$seed" "$job_label" "$state" "$elapsed_h" "${pct}%" "${current_ns}ns" "$ns_day" "$eta_h" "$slow_flag"
done

echo
echo "Target: 200 ns (step $TARGET_STEPS). Benchmark throughput (with --exclusive,"
echo "the original single-system run): ${BENCHMARK_NS_DAY} ns/day. Since these jobs"
echo "run WITHOUT --exclusive, ns/day well below that benchmark is expected to be"
echo "worth investigating, not necessarily alarming on its own -- but if it's"
echo "dramatically lower (e.g. <10% of benchmark), check for node contention."

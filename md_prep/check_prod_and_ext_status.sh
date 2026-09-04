#!/bin/bash
#
# Machine-readable status for the prod/extension job pair of each of the 8
# seeds. Companion to check_prod_performance.sh (that one is for human-
# readable throughput monitoring; this one is for scripted polling that
# needs to detect state transitions -- prod job TIMEOUT + extension job
# START, and target reached / 200 ns complete).
#
# Run this ON THE CLUSTER (uses sacct and /scratch/alpine paths).
#
# Prints one pipe-delimited line per seed:
#   seed|prod_state|prod_exit|ext_state|ext_exit|last_step|current_ns|pct|target_reached
#
# target_reached is "yes" if last_step >= TARGET_STEPS, else "no".
# Any field is "-" if not yet available (e.g. ext job hasn't started).
#
# Usage: bash check_prod_and_ext_status.sh [seed ...]

SCRATCH_ROOT=/scratch/alpine/ivta1597/SH2/seeds
ALL_SEEDS="af3_c1 af3_c2 af3_c3 af3_c4 boltz_c1 boltz_c2 boltz_c3 boltz_c4"
SEEDS="${*:-$ALL_SEEDS}"

TARGET_STEPS=100000000
DT_NS=0.000002

sacct_job () {
  # $1 = exact job name -- prints "State Exit JobID" for the most recent
  # matching job (tab-separated), or empty if none found
  sacct -u "$USER" --name="$1" --format=State,ExitCode,JobID --noheader -X 2>/dev/null | tail -1
}

squeue_job () {
  # $1 = exact job name -- prints "State -" (sacct-compatible, exit code
  # unknown while pending/running) for a job still in the queue (covers
  # PENDING jobs held on a dependency, which sacct doesn't show until they
  # actually start)
  local st
  st=$(squeue -u "$USER" --name="$1" --noheader --format="%T" 2>/dev/null | head -1)
  [ -n "$st" ] && echo "$st -"
}

for seed in $SEEDS; do
  prod_info=$(sacct_job "${seed}_prod")
  ext_info=$(sacct_job "${seed}_prod_xtnd")
  [ -z "$ext_info" ] && ext_info=$(squeue_job "${seed}_prod_xtnd")

  prod_state=$(echo "$prod_info" | awk '{print $1}'); prod_state="${prod_state:--}"
  prod_exit=$(echo "$prod_info" | awk '{print $2}'); prod_exit="${prod_exit:--}"
  ext_state=$(echo "$ext_info" | awk '{print $1}'); ext_state="${ext_state:--}"
  ext_exit=$(echo "$ext_info" | awk '{print $2}'); ext_exit="${ext_exit:--}"

  DIR="$SCRATCH_ROOT/$seed/prod_md_48x1"
  log="$DIR/prod_md_200ns.log"

  last_step="-"; current_ns="-"; pct="-"; target_reached="no"
  if [ -f "$log" ]; then
    last_step=$(grep -oE "Writing checkpoint, step [0-9]+" "$log" | tail -1 | awk '{print $NF}')
    if [ -z "$last_step" ]; then
      last_step=$(grep -A1 "^ *Step *Time" "$log" | tail -1 | awk '{print $1}')
    fi
    if [ -n "$last_step" ]; then
      current_ns=$(echo "scale=3; $last_step * $DT_NS" | bc)
      pct=$(echo "scale=1; $last_step * 100 / $TARGET_STEPS" | bc)
      if [ "$last_step" -ge "$TARGET_STEPS" ]; then
        target_reached="yes"
      fi
    else
      last_step="-"
    fi
  fi

  echo "${seed}|${prod_state}|${prod_exit}|${ext_state}|${ext_exit}|${last_step}|${current_ns}|${pct}|${target_reached}"
done

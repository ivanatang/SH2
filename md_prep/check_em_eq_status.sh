#!/bin/bash
#
# Check that EM and equilibration (NVT+NPT) finished as expected for all 8
# seeds. Run this ON THE CLUSTER (uses /scratch/alpine paths, gmx, and sacct)
# after submitting em_<seed>.sh / equil_<seed>.sh.
#
# For each seed, per stage, checks:
#   - SLURM job state via sacct (by job name), if any job with that name exists
#   - whether the expected final .gro was written (mdrun only writes -c output
#     on clean completion -- a killed/crashed run leaves a .cpt but no .gro)
#   - "Fatal error" in the GROMACS log
#   - EM: converged to Fmax < 1000, or ran to max steps without converging
#   - NVT: average temperature (expect ~300 K)
#   - NPT: average pressure (expect ~1 bar, noisy over only 100 ps) and density
#
# Usage: bash check_em_eq_status.sh [seed ...]
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

avg_energy_term () {
  # $1 = .edr file, $2 = term name (e.g. "Temperature", "Pressure", "Density")
  local edr="$1" term="$2"
  [ -f "$edr" ] || { echo "N/A"; return; }
  echo "$term
0" | gmx energy -f "$edr" -o /tmp/_check_em_eq_$$.xvg 2>/dev/null \
    | grep -i "^${term}" | head -1 | awk '{print $2}'
  rm -f /tmp/_check_em_eq_$$.xvg
}

sacct_state () {
  # $1 = job name -- most recent matching job's state, or "no job found"
  local name="$1"
  local line
  line=$(sacct -u "$USER" --name="$name" --format=JobID,State --noheader -X 2>/dev/null | tail -1)
  if [ -z "$line" ]; then
    echo "no job found"
  else
    echo "$line" | awk '{print $2}'
  fi
}

printf "%-10s | %-8s %-6s %-10s %-12s | %-8s %-6s %-10s | %-8s %-6s %-10s %-10s\n" \
  "seed" "EM:sacct" "gro?" "Fmax" "converged?" "NVT:sacct" "gro?" "avgT(K)" "NPT:sacct" "gro?" "avgP(bar)" "density"
printf '%.0s-' {1..140}; echo

for seed in $SEEDS; do
  DIR="$SCRATCH_ROOT/$seed"

  # ---------------- EM ----------------
  em_sacct=$(sacct_state "em_${seed}")
  em_gro="no"; [ -f "$DIR/EM/em.gro" ] && em_gro="yes"
  em_fmax="N/A"; em_conv="N/A"
  if [ -f "$DIR/EM/em.log" ]; then
    if grep -qi "Fatal error" "$DIR/EM/em.log"; then
      em_conv="FATAL_ERROR"
    elif grep -qi "converged to Fmax" "$DIR/EM/em.log"; then
      em_conv="yes"
    elif grep -qi "did not converge" "$DIR/EM/em.log"; then
      em_conv="no(maxsteps)"
    fi
    em_fmax=$(grep -i "Maximum force" "$DIR/EM/em.log" | tail -1 | awk '{print $4}')
  fi

  # ---------------- NVT ----------------
  nvt_sacct=$(sacct_state "eq_${seed}")
  nvt_gro="no"; [ -f "$DIR/NVT/nvt.gro" ] && nvt_gro="yes"
  nvt_err="ok"
  [ -f "$DIR/NVT/nvt.log" ] && grep -qi "Fatal error" "$DIR/NVT/nvt.log" && nvt_err="FATAL_ERROR"
  nvt_temp="N/A"
  [ "$nvt_err" = "ok" ] && nvt_temp=$(avg_energy_term "$DIR/NVT/nvt.edr" "Temperature")

  # ---------------- NPT ----------------
  npt_sacct="$nvt_sacct"  # same job (equil_<seed>.sh runs NVT then NPT)
  npt_gro="no"; [ -f "$DIR/NPT/npt.gro" ] && npt_gro="yes"
  npt_err="ok"
  [ -f "$DIR/NPT/npt.log" ] && grep -qi "Fatal error" "$DIR/NPT/npt.log" && npt_err="FATAL_ERROR"
  npt_pres="N/A"; npt_dens="N/A"
  if [ "$npt_err" = "ok" ]; then
    npt_pres=$(avg_energy_term "$DIR/NPT/npt.edr" "Pressure")
    npt_dens=$(avg_energy_term "$DIR/NPT/npt.edr" "Density")
  fi

  printf "%-10s | %-8s %-6s %-10s %-12s | %-8s %-6s %-10s | %-8s %-6s %-10s %-10s\n" \
    "$seed" "$em_sacct" "$em_gro" "$em_fmax" "$em_conv" \
    "$nvt_sacct" "$nvt_gro" "$nvt_temp" \
    "$npt_sacct" "$npt_gro" "$npt_pres" "$npt_dens"
done

echo
echo "Expect: EM converged=yes (or no(maxsteps) is non-fatal but worth a look),"
echo "        NVT avgT ~300 K, NPT avgP ~1 bar (noisy over just 100 ps) and a"
echo "        stable density. Anything FATAL_ERROR or gro=no needs investigating"
echo "        before launching production from that seed."

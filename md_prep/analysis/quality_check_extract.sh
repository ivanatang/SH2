#!/bin/bash
#SBATCH --job-name=qc_extract
#SBATCH --output=qc_extract_%j.out
#SBATCH --error=qc_extract_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --qos=cpu-normal
#
# Stage 1 of the production-trajectory quality check (companion: quality_check_assess.py).
# For each of the 8 seeds, extracts the raw signals needed to judge:
#   - atom clashes / blow-ups   (log warnings + potential energy trace)
#   - RMSD trace shape          (protein backbone RMSD to the run's start structure)
#   - peptide dissociation      (min distance + contact count, protein vs peptide)
#
# Does NOT judge pass/fail itself -- quality_check_assess.py reads what this writes
# and applies the actual thresholds, so the numeric criteria live in one place.
#
# Group names ("Protein", "LIG", "Backbone") were confirmed identical across all
# 8 seeds via `gmx_mpi make_ndx` on 2026-09-04 -- Protein is the 112-residue protein
# chain (1778 atoms), LIG is the pTyr peptide, parameterized as a single ligand-like
# residue rather than a standard polymer (hence it doesn't show up as "Protein").
#
# Usage: sbatch quality_check_extract.sh [seed ...]
#   (no args = all 8 seeds; or pass specific seeds, e.g. af3_c1 boltz_c2)

export PATH=/curc/sw/install/openmpi/5.0.6/gcc/14.2.0/bin:$PATH
export LD_LIBRARY_PATH=/curc/sw/install/gcc/14.2.0/lib64:/curc/sw/install/openmpi/5.0.6/gcc/14.2.0/lib:$LD_LIBRARY_PATH
source /curc/sw/install/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate SH2

module purge
module load gromacs
# gmx_mpi (not plain gmx) is required -- see md_prep/seeds/*/resume_prod_*.sh comments
# for why the gcc/openmpi/anaconda modules are bypassed here.

SCRATCH_ROOT=/scratch/alpine/ivta1597/SH2/seeds
ALL_SEEDS="af3_c1 af3_c2 af3_c3 af3_c4 boltz_c1 boltz_c2 boltz_c3 boltz_c4"
SEEDS="${*:-$ALL_SEEDS}"

CONTACT_CUTOFF_NM=0.5   # heavy-atom contact distance for the peptide dissociation check

for seed in $SEEDS; do
  DIR="$SCRATCH_ROOT/$seed/prod_md_48x1"
  QC="$DIR/qc"
  mkdir -p "$QC"

  tpr="$DIR/prod_md_200ns.tpr"
  xtc="$DIR/prod_md_200ns.xtc"
  edr="$DIR/prod_md_200ns.edr"
  log="$DIR/prod_md_200ns.log"

  if [ ! -f "$tpr" ] || [ ! -f "$xtc" ]; then
    echo "$seed: missing $tpr or $xtc -- skipping"
    continue
  fi

  echo "=== $seed ==="

  # 1. Log warning signatures for clashes/blow-ups/instability.
  if [ ! -f "$QC/log_warnings_count.txt" ]; then
    grep -nE "Fatal error|LINCS WARNING|1-4 interaction|can not be settled|Water molecule starting|pressure scaling|maximum distance between bonded" \
      "$log" > "$QC/log_warnings.txt" 2>/dev/null
    wc -l < "$QC/log_warnings.txt" > "$QC/log_warnings_count.txt"
  fi

  # 2. Potential/total energy + temperature trace (NaN/Inf or spikes -> blow-up).
  if [ ! -s "$QC/energy.xvg" ]; then
    printf "Potential\nTotal-Energy\nTemperature\n\n" | \
      gmx_mpi energy -s "$tpr" -f "$edr" -o "$QC/energy.xvg" -xvg none
  fi

  # 3. Protein backbone RMSD to the tpr's own reference (t=0 of this 200ns target).
  if [ ! -s "$QC/rmsd_backbone.xvg" ]; then
    printf "Backbone\nBackbone\n" | \
      gmx_mpi rms -s "$tpr" -f "$xtc" -o "$QC/rmsd_backbone.xvg" -tu ns -xvg none
  fi

  # 4. Protein-peptide minimum distance + contact count (dissociation check).
  if [ ! -s "$QC/mindist_prot_lig.xvg" ] || [ ! -s "$QC/numcontacts_prot_lig.xvg" ]; then
    printf "Protein\nLIG\n" | \
      gmx_mpi mindist -s "$tpr" -f "$xtc" -od "$QC/mindist_prot_lig.xvg" \
        -on "$QC/numcontacts_prot_lig.xvg" -d "$CONTACT_CUTOFF_NM" -tu ns -xvg none
  fi

  echo "$seed: done -> $QC"
done

echo
echo "Extraction complete. Run quality_check_assess.py against these seed dirs next."

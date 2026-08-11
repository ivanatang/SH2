#!/bin/bash
# postprocess_SH2.sh
# ─────────────────────────────────────────────────────────────────────────────
# Post-processes the SH2 superbinder + pTyr peptide production trajectory:
# PBC correction, fitting, RMSD/RMSF-style QC, and the project-specific
# "did the phosphate stay in the engineered pocket" contact metric (the SH2
# analog of biosensors' gate-latch distance).
#
# This runs locally against the OneDrive-synced run directory (no SLURM/Alpine
# needed for a single ~50 ns, ~21k-atom trajectory) -- unlike biosensors'
# multi-sequence post_processing_pipeline_worker.sh, there's no config.yaml or
# per-sequence looping here, just one system.
#
# Pocket residues (40, 45, 63 in this system's numbering) were verified via
# md_prep/gromacs/sh2_superbinder_pTyr_dodecahedron_protein.itp to be
# VAL/ALA/LEU -- i.e. the T183V/C188A/K206L engineered superbinder mutations
# (UniProt_resi - 143 = this system's resi, same offset established in
# af3/README.md and used in af3/align_af3_vs_crystal.py's MUTATED_POSITIONS).
#
# Usage:
#   bash postprocess_SH2.sh [DATA_DIR] [FORCE]
#     DATA_DIR  directory containing EM/ NVT/ NPT/ prod_md_48x1/ subdirs
#               (default: the OneDrive-synced SH2_enhanced_sampling folder)
#     FORCE     true to overwrite existing outputs (default: false)
#
# Outputs land in $DATA_DIR/analysis/:
#   energy_*.xvg                        EM/NVT/NPT/production energy checks
#   index.ndx                           custom groups (Protein_Peptide, Peptide_heavy,
#                                        Phosphate, Pocket_muts_heavy)
#   pbc.xtc, fit.xtc                    PBC-corrected / fitted, whole system
#   PL_only.xtc                         protein+peptide only (solvent stripped),
#                                        used for the movie and lighter analyses
#   rmsd_protein_backbone.xvg           backbone RMSD vs t=0 (equilibration/stability QC)
#   rmsd_peptide_heavy.xvg              peptide heavy-atom RMSD vs t=0 (pose stability)
#   peptide_protein_mindist.xvg         min distance peptide<->protein (dissociation check)
#   phosphate_pocket_mindist.xvg        min distance phosphate<->pocket mutations (KEY metric)
#   phosphate_pocket_comdist.xvg        COM-COM distance, same pair (smoother companion signal)
#   Rg_complex.xvg                      radius of gyration of the complex (global compactness)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

GMX=${GMX:-gmx}

DEFAULT_DATA_DIR="/Users/ivanatang/Library/CloudStorage/OneDrive-UCB-O365/Shirts Lab/SH2_enhanced_sampling"
DATA_DIR="${1:-$DEFAULT_DATA_DIR}"
FORCE="${2:-false}"

if [[ "$FORCE" != "true" && "$FORCE" != "false" ]]; then
    echo "ERROR: FORCE must be true or false (got: '$FORCE')" >&2; exit 1
fi
should_run() { [[ "$FORCE" == "true" || ! -f "$1" ]]; }

EM_DIR="$DATA_DIR/EM"
NVT_DIR="$DATA_DIR/NVT"
NPT_DIR="$DATA_DIR/NPT"
PROD_DIR="$DATA_DIR/prod_md_48x1"
TPR="$PROD_DIR/prod_md_500ns.tpr"
XTC="$PROD_DIR/prod_md_500ns.xtc"
OUT_DIR="$DATA_DIR/analysis"

for f in "$TPR" "$XTC"; do
    [[ -f "$f" ]] || { echo "ERROR: required file not found: $f" >&2; exit 1; }
done

mkdir -p "$OUT_DIR"
NDX="$OUT_DIR/index.ndx"
PBC_XTC="$OUT_DIR/pbc.xtc"
FIT_XTC="$OUT_DIR/fit.xtc"
PL_XTC="$OUT_DIR/PL_only.xtc"

echo "════════════════════════════════════════════════════════════════"
echo "  SH2 postprocessing: $DATA_DIR"
echo "  Output directory: $OUT_DIR"
echo "════════════════════════════════════════════════════════════════"

# ── Step 0: Energy checks (EM/NVT/NPT/production) ─────────────────────────────
# Term names (not numbers) are used deliberately -- gmx energy matches them
# directly, and numeric term indices shift depending on which terms a given
# .edr happens to contain (confirmed different between EM/NVT/NPT/production
# for this system), so hardcoding numbers the way biosensors' run_energy_check.sh
# does would be fragile here.
echo "[0/11] Energy checks"
if should_run "$OUT_DIR/energy_em_potential.xvg"; then
    echo "  EM potential → energy_em_potential.xvg"
    echo "Potential" | "$GMX" energy -f "$EM_DIR/em.edr" -o "$OUT_DIR/energy_em_potential.xvg" > "$OUT_DIR/energy_em_potential.log" 2>&1
fi
if should_run "$OUT_DIR/energy_nvt_temperature.xvg"; then
    echo "  NVT temperature → energy_nvt_temperature.xvg"
    echo "Temperature" | "$GMX" energy -f "$NVT_DIR/nvt.edr" -o "$OUT_DIR/energy_nvt_temperature.xvg" > "$OUT_DIR/energy_nvt_temperature.log" 2>&1
fi
if should_run "$OUT_DIR/energy_npt_pressure.xvg"; then
    echo "  NPT pressure → energy_npt_pressure.xvg"
    echo "Pressure" | "$GMX" energy -f "$NPT_DIR/npt.edr" -o "$OUT_DIR/energy_npt_pressure.xvg" > "$OUT_DIR/energy_npt_pressure.log" 2>&1
fi
if should_run "$OUT_DIR/energy_npt_density.xvg"; then
    echo "  NPT density → energy_npt_density.xvg"
    echo "Density" | "$GMX" energy -f "$NPT_DIR/npt.edr" -o "$OUT_DIR/energy_npt_density.xvg" > "$OUT_DIR/energy_npt_density.log" 2>&1
fi
if should_run "$OUT_DIR/energy_prod_potential.xvg"; then
    echo "  Production potential → energy_prod_potential.xvg"
    echo "Potential" | "$GMX" energy -f "$PROD_DIR/prod_md_500ns.edr" -o "$OUT_DIR/energy_prod_potential.xvg" > "$OUT_DIR/energy_prod_potential.log" 2>&1
fi
if should_run "$OUT_DIR/energy_prod_temperature.xvg"; then
    echo "  Production temperature → energy_prod_temperature.xvg"
    echo "Temperature" | "$GMX" energy -f "$PROD_DIR/prod_md_500ns.edr" -o "$OUT_DIR/energy_prod_temperature.xvg" > "$OUT_DIR/energy_prod_temperature.log" 2>&1
fi
if should_run "$OUT_DIR/energy_prod_pressure.xvg"; then
    echo "  Production pressure → energy_prod_pressure.xvg"
    echo "Pressure" | "$GMX" energy -f "$PROD_DIR/prod_md_500ns.edr" -o "$OUT_DIR/energy_prod_pressure.xvg" > "$OUT_DIR/energy_prod_pressure.log" 2>&1
fi
if should_run "$OUT_DIR/energy_prod_density.xvg"; then
    echo "  Production density → energy_prod_density.xvg"
    echo "Density" | "$GMX" energy -f "$PROD_DIR/prod_md_500ns.edr" -o "$OUT_DIR/energy_prod_density.xvg" > "$OUT_DIR/energy_prod_density.log" 2>&1
fi

# ── Step 1: Custom index groups ────────────────────────────────────────────────
# Peptide_heavy      peptide (resname LIG) minus hydrogens
# Protein_Peptide    protein + peptide, used for centering/fitting
# Phosphate          the pTyr phosphate group (P, O1P, O2P, O3P)
# Pocket_muts_heavy  heavy atoms of residues 40/45/63 (= V183/A188/L206 superbinder muts)
if should_run "$NDX"; then
    echo "[1/11] make_ndx → $NDX"
    # Staged (no parentheses): gmx make_ndx's expression parser rejects
    # "(ri 40 | ri 45 | ri 63) & ! a H*" as a single expression -- confirmed
    # by testing directly against this system's .tpr. Build the union group
    # first (index 23, auto-numbered after the 3 groups above it), then AND
    # it against "not hydrogen" as a separate step (index 24) referencing it
    # by its now-deterministic numeric index (named references to a
    # just-renamed group also failed parsing; numeric references work).
    "$GMX" make_ndx -f "$TPR" -o "$NDX" << 'EOF'
r LIG & ! a H*
name 20 Peptide_heavy
1 | r LIG
name 21 Protein_Peptide
r LIG & a P O1P O2P O3P
name 22 Phosphate
ri 40 | ri 45 | ri 63
name 23 Pocket_muts
23 & ! a H*
name 24 Pocket_muts_heavy
q
EOF
else
    echo "[1/11] SKIP make_ndx ($NDX exists)"
fi

# ── Step 2: PBC correction ─────────────────────────────────────────────────────
if should_run "$PBC_XTC"; then
    echo "[2/11] trjconv PBC → $PBC_XTC"
    echo "Protein_Peptide System" | "$GMX" trjconv \
        -s "$TPR" -f "$XTC" -n "$NDX" -o "$PBC_XTC" \
        -pbc mol -center
else
    echo "[2/11] SKIP PBC trjconv ($PBC_XTC exists)"
fi

# ── Step 3: Rotational/translational fitting ───────────────────────────────────
if should_run "$FIT_XTC"; then
    echo "[3/11] trjconv fit → $FIT_XTC"
    echo "Protein_Peptide System" | "$GMX" trjconv \
        -s "$TPR" -f "$PBC_XTC" -n "$NDX" -o "$FIT_XTC" \
        -fit rot+trans
else
    echo "[3/11] SKIP fit trjconv ($FIT_XTC exists)"
fi

# ── Step 4: Extract protein+peptide-only trajectory (solvent stripped) ────────
if should_run "$PL_XTC"; then
    echo "[4/11] Extract Protein_Peptide-only → $PL_XTC"
    echo "Protein_Peptide" | "$GMX" trjconv \
        -s "$TPR" -f "$FIT_XTC" -n "$NDX" -o "$PL_XTC"
else
    echo "[4/11] SKIP PL_only ($PL_XTC exists)"
fi

# ── Step 5: Protein backbone RMSD vs t=0 ───────────────────────────────────────
# Equilibration/stability QC: did the protein fold stay put over 50 ns?
if should_run "$OUT_DIR/rmsd_protein_backbone.xvg"; then
    echo "[5/11] rms (Backbone vs t=0) → rmsd_protein_backbone.xvg"
    echo "Backbone Backbone" | "$GMX" rms \
        -s "$TPR" -f "$FIT_XTC" -n "$NDX" \
        -o "$OUT_DIR/rmsd_protein_backbone.xvg"
else
    echo "[5/11] SKIP rmsd_protein_backbone.xvg (exists)"
fi

# ── Step 6: Peptide heavy-atom RMSD vs t=0 ─────────────────────────────────────
# Pose stability: is the peptide holding a consistent conformation/pose, or
# drifting/refolding within the pocket?
if should_run "$OUT_DIR/rmsd_peptide_heavy.xvg"; then
    echo "[6/11] rms (Peptide_heavy vs t=0) → rmsd_peptide_heavy.xvg"
    echo "Peptide_heavy Peptide_heavy" | "$GMX" rms \
        -s "$TPR" -f "$FIT_XTC" -n "$NDX" \
        -o "$OUT_DIR/rmsd_peptide_heavy.xvg"
else
    echo "[6/11] SKIP rmsd_peptide_heavy.xvg (exists)"
fi

# ── Step 7: Peptide<->protein minimum distance ─────────────────────────────────
# Dissociation check: does this ever depart from van-der-Waals contact
# (~0.3-0.4 nm), which would indicate the peptide unbound?
#
# IMPORTANT: computed on the RAW trajectory ($XTC), not $FIT_XTC. gmx trjconv
# -fit rot+trans rotates coordinates but does NOT rotate the stored box
# vectors to match. gmx mindist/gmx distance apply their own periodic-image
# correction using those box vectors by default, so on a rotated-but-not-
# rebox'd trajectory they can spuriously wrap an atom into the wrong image,
# producing large fake "jumps" that look like transient dissociation but
# aren't real -- confirmed by direct comparison: the same window computed on
# $FIT_XTC showed wild 0.02-0.65 nm swings where the raw trajectory shows a
# flat 0.25-0.28 nm. The raw trajectory's box vectors are always internally
# consistent with its coordinates, so mindist's built-in PBC handling (-pbc
# yes, default) is correct there without needing pbc.xtc/fit.xtc at all.
if should_run "$OUT_DIR/peptide_protein_mindist.xvg"; then
    echo "[7/11] mindist Peptide_heavy<->Protein-H (raw trajectory) → peptide_protein_mindist.xvg"
    echo "Peptide_heavy Protein-H" | "$GMX" mindist \
        -s "$TPR" -f "$XTC" -n "$NDX" \
        -od "$OUT_DIR/peptide_protein_mindist.xvg"
else
    echo "[7/11] SKIP peptide_protein_mindist.xvg (exists)"
fi

# ── Step 8: Phosphate<->pocket-mutation minimum distance (KEY METRIC) ─────────
# This is the SH2-project analog of biosensors' gate-latch distance: the
# primary observable for whether the peptide's phosphate stayed engaged in
# the engineered superbinder pocket (V183/A188/L206) throughout production,
# vs. drifting away or dissociating.
#
# Computed on the RAW trajectory -- see step 7's comment for why.
if should_run "$OUT_DIR/phosphate_pocket_mindist.xvg"; then
    echo "[8/11] mindist Phosphate<->Pocket_muts_heavy (raw trajectory) → phosphate_pocket_mindist.xvg"
    echo "Phosphate Pocket_muts_heavy" | "$GMX" mindist \
        -s "$TPR" -f "$XTC" -n "$NDX" \
        -od "$OUT_DIR/phosphate_pocket_mindist.xvg"
else
    echo "[8/11] SKIP phosphate_pocket_mindist.xvg (exists)"
fi

# ── Step 9: Phosphate<->pocket-mutation center-of-mass distance ───────────────
# Smoother companion signal to the mindist above (mindist can be noisy frame
# to frame since it tracks whichever single atom pair is closest at each step).
#
# Computed on the RAW trajectory -- see step 7's comment for why.
if should_run "$OUT_DIR/phosphate_pocket_comdist.xvg"; then
    echo "[9/11] distance (COM) Phosphate<->Pocket_muts_heavy (raw trajectory) → phosphate_pocket_comdist.xvg"
    "$GMX" distance \
        -s "$TPR" -f "$XTC" -n "$NDX" \
        -select 'com of group "Phosphate" plus com of group "Pocket_muts_heavy"' \
        -oall "$OUT_DIR/phosphate_pocket_comdist.xvg"
else
    echo "[9/11] SKIP phosphate_pocket_comdist.xvg (exists)"
fi

# ── Step 10: Radius of gyration of the complex ─────────────────────────────────
# Sanity check for gross unfolding/instability over the trajectory.
if should_run "$OUT_DIR/Rg_complex.xvg"; then
    echo "[10/11] gyrate (Protein_Peptide) → Rg_complex.xvg"
    echo "Protein_Peptide" | "$GMX" gyrate \
        -s "$TPR" -f "$FIT_XTC" -n "$NDX" \
        -o "$OUT_DIR/Rg_complex.xvg"
else
    echo "[10/11] SKIP Rg_complex.xvg (exists)"
fi

# ── Step 11: Clean up GROMACS backup files ─────────────────────────────────────
echo "[11/11] Removing GROMACS backup files (#*#)"
rm -f "$DATA_DIR"/\#*\# "$OUT_DIR"/\#*\#

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  Postprocessing complete: $OUT_DIR"
echo "════════════════════════════════════════════════════════════════"

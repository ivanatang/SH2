"""
align_af3_vs_crystal.py

Structural validation: does AlphaFold3's predicted phosphopeptide pose engage
the REAL superbinder SH2 pocket seen in the 4F5B crystal structure, or did it
just dock the peptide somewhere generically plausible on the protein surface?

Inputs
------
1. Crystal structure (ground truth pocket definition):
     4F5B - triple-mutant (T183V/C188A/K206L) c-Src SH2 domain bound to FREE
     phosphotyrosine. Chain A = protein, resn PTR / resi 301 / chain A = ligand.

2. AlphaFold3 model (prediction to be validated):
     Same SH2 domain (chain A) predicted in complex with an 11-mer
     phosphopeptide DHEPIYEQWGW (chain B), where residue 6 (Tyr) was specified
     as phosphotyrosine (resn PTR, chain B, resi 6).

Logic
-----
A. Superpose AF3 chain A onto crystal chain A (protein only -- same protein
   sequence in both, so this is a same-sequence structural sanity check, not
   a homology alignment). PyMOL's cmd.align() transforms the ENTIRE mobile
   object (not just the aligned atoms), so this rigid-body transform carries
   the AF3-predicted peptide (chain B) along with it into the crystal's
   reference frame. A low RMSD here (~<2 A for a well-resolved,
   well-characterized domain like this) is just a sanity check that AF3 got
   the fold right -- it says nothing about whether the peptide is docked
   correctly.

B. Once both structures share a reference frame, we compare where the
   AF3-predicted phosphate group landed vs. where the crystallographic free
   phosphotyrosine's phosphate actually sits. The phosphorus atom (atom name
   "P") is used as the comparison point because it is the single atom most
   diagnostic of pocket engagement -- the pTyr phosphate is what is
   recognized by the SH2 domain's conserved Arg-beta5 (the phosphate-binding
   pocket), so if AF3 got the phosphate location right, it almost certainly
   got the key specificity contact right too.

C. We independently define "the real pocket" as every protein residue within
   POCKET_CUTOFF (default 5 A) of the crystallographic PTR 301 -- this should
   include the three engineered superbinder mutations (183/188/206) plus
   other native contacts (e.g. the conserved Arg/Ser of the phosphate pocket).
   We then check, in the SAME reference frame, whether the AF3-predicted PTR
   is also within POCKET_CUTOFF of those SAME residues. This is a stronger
   test than the raw P-P distance alone: a molecule could sit close in space
   without touching the actual contact-forming residues (e.g. it could be
   float above the pocket mouth), so both checks are reported.

Robustness notes (why this is a .py script, not .pml)
------------------------------------------------------
AlphaFold3 mmCIF output uses standard mmCIF label_*/auth_* dual numbering.
PyMOL's cif loader uses auth_asym_id for `chain` and auth_seq_id for `resi`,
which for this AF3 job happens to line up with the PDB-style chain B / resi 6
convention requested. However, mmCIF outputs from different tools/versions
can vary (some use "P1" instead of "P" for the phosphate phosphorus, chain
letters can shift, etc.), so every selection that MUST find exactly one atom
is wrapped in `find_atom()` below, which tries a list of candidate atom names
and fails loudly with a diagnostic atom dump (rather than silently operating
on an empty selection) if none match. If you see a "COULD NOT FIND ATOM"
error, open the .cif in a text editor, find the PTR HETATM block, and add
the actual atom name you see to the relevant *_P_ATOM_CANDIDATES list below.

Usage
-----
    pymol -cq /Users/ivanatang/Developer/SH2/af3/align_af3_vs_crystal.py

Before running, confirm the AF3 model filename (AF3 sometimes nests the
representative model under a job-name prefix, or you may want a specific
seed/sample instead of the top-ranked model). To check:
    ls -la /Users/ivanatang/Developer/SH2/af3/output/superbinder_sh2_ptyr_peptide/
The script auto-detects a reasonable file (see resolve_af3_cif() below) and
prints exactly which file it used -- check the log for that line.

Outputs
-------
    <OUTDIR>/af3_vs_crystal_validation.pse   (session, both structures aligned)
    <OUTDIR>/af3_vs_crystal_pocket.png       (rendered close-up of the pocket)
plus all numeric results printed to stdout/log.
"""

import glob
import os
import sys

from pymol import cmd

# --------------------------------------------------------------------------
# CONFIG -- edit these if your paths, chain IDs, or residue numbers differ.
# --------------------------------------------------------------------------
CRYSTAL_PATH = "/Users/ivanatang/Developer/SH2/structures/4F5B_superbinder_SH2_pTyr.pdb"
AF3_DIR = "/Users/ivanatang/Developer/SH2/af3/output/superbinder_sh2_ptyr_peptide"
AF3_PREFERRED_NAME = "superbinder_sh2_ptyr_peptide_model.cif"  # AF3's top-ranked model
OUTDIR = "/Users/ivanatang/Developer/SH2/af3"

CRYSTAL_PROT_CHAIN = "A"
CRYSTAL_LIG_CHAIN = "A"
CRYSTAL_LIG_RESI = "301"

AF3_PROT_CHAIN = "A"
AF3_PEPTIDE_CHAIN = "B"
AF3_LIG_RESI = "6"

LIG_RESN = "PTR"
P_ATOM_CANDIDATES = ["P", "P1"]        # phosphate phosphorus, name varies by source
MUTATED_POSITIONS = ["183", "188", "206"]  # T183V / C188A / K206L (superbinder muts)
POCKET_CUTOFF = 6.0                    # Angstroms, defines "real pocket" contacts
REMOVE_SOLVENT = True                  # strip crystallographic waters for clarity

SESSION_OUT = os.path.join(OUTDIR, "af3_vs_crystal_validation.pse")
PNG_OUT = os.path.join(OUTDIR, "af3_vs_crystal_pocket.png")


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def resolve_af3_cif():
    """Find the AF3 model .cif, tolerating AF3's variable output naming."""
    preferred = os.path.join(AF3_DIR, AF3_PREFERRED_NAME)
    if os.path.isfile(preferred):
        print(f"[INFO] Using AF3 model: {preferred}")
        return preferred

    print(f"[WARN] Expected AF3 file not found at {preferred}")
    print(f"[WARN] Searching {AF3_DIR} for alternatives...")

    candidates = sorted(glob.glob(os.path.join(AF3_DIR, "*_model.cif")))
    if not candidates:
        candidates = sorted(glob.glob(os.path.join(AF3_DIR, "**", "model.cif"), recursive=True))
    if not candidates:
        candidates = sorted(glob.glob(os.path.join(AF3_DIR, "**", "*.cif"), recursive=True))

    if not candidates:
        sys.exit(
            f"[FATAL] No .cif files found under {AF3_DIR}. "
            f"Confirm AF3 output has been downloaded and check AF3_DIR at the "
            f"top of this script."
        )

    chosen = candidates[0]
    print(f"[INFO] Falling back to: {chosen}")
    print(f"[INFO] (All candidates found: {candidates})")
    return chosen


def find_atom(base_selection, atom_candidates, label):
    """
    Return a selection string that resolves to exactly one atom, trying each
    name in atom_candidates in turn. Exits with a diagnostic atom dump if
    none match (rather than silently returning an empty/ambiguous selection).
    """
    for name in atom_candidates:
        sel = f"({base_selection}) and name {name}"
        n = cmd.count_atoms(sel)
        if n == 1:
            print(f"[INFO] {label}: found atom '{name}' via selection: {sel}")
            return sel
        if n > 1:
            print(f"[WARN] {label}: name '{name}' matched {n} atoms (expected 1), skipping.")

    # Nothing matched -- dump what IS present so the user can fix the config.
    print(f"[ERROR] {label}: could not find a unique atom among {atom_candidates}.")
    print(f"[ERROR] Atoms present in base selection '{base_selection}':")
    atom_names = set()
    cmd.iterate(base_selection, "atom_names.add(name)", space={"atom_names": atom_names})
    if atom_names:
        print(f"        {sorted(atom_names)}")
    else:
        print("        (base selection is EMPTY -- check chain/resi/resn in CONFIG above)")
    sys.exit(
        f"[FATAL] Add the correct atom name to the candidate list for {label} "
        f"and re-run."
    )


def require_atoms(selection, label):
    n = cmd.count_atoms(selection)
    if n == 0:
        sys.exit(
            f"[FATAL] Selection for {label} is EMPTY: '{selection}'. "
            f"Check chain IDs / resi numbers in the CONFIG block."
        )
    return n


# --------------------------------------------------------------------------
# 1. Load structures
# --------------------------------------------------------------------------
if not os.path.isfile(CRYSTAL_PATH):
    sys.exit(f"[FATAL] Crystal structure not found: {CRYSTAL_PATH}")

af3_path = resolve_af3_cif()

cmd.load(CRYSTAL_PATH, "crystal")
cmd.load(af3_path, "af3_model")

if REMOVE_SOLVENT:
    cmd.remove("crystal and solvent")
    cmd.remove("af3_model and solvent")

require_atoms(f"crystal and chain {CRYSTAL_PROT_CHAIN} and polymer.protein", "crystal protein chain A")
require_atoms(f"af3_model and chain {AF3_PROT_CHAIN} and polymer.protein", "AF3 protein chain A")
require_atoms(
    f"crystal and chain {CRYSTAL_LIG_CHAIN} and resi {CRYSTAL_LIG_RESI} and resn {LIG_RESN}",
    "crystal free pTyr ligand (PTR 301)",
)
require_atoms(
    f"af3_model and chain {AF3_PEPTIDE_CHAIN} and resi {AF3_LIG_RESI} and resn {LIG_RESN}",
    "AF3-predicted peptide pTyr (chain B resi 6)",
)

# --------------------------------------------------------------------------
# 2. Superpose AF3 chain A onto crystal chain A (protein backbone/sequence
#    alignment). This moves the ENTIRE af3_model object -- including the
#    predicted chain B peptide -- into the crystal's reference frame.
# --------------------------------------------------------------------------
mobile_sel = f"af3_model and chain {AF3_PROT_CHAIN} and polymer.protein"
target_sel = f"crystal and chain {CRYSTAL_PROT_CHAIN} and polymer.protein"

align_result = cmd.align(mobile_sel, target_sel, object="aln_chainA")
rmsd_after, n_atoms_after, n_cycles, rmsd_before, n_atoms_before, score, n_res_aligned = align_result

print("\n" + "=" * 70)
print("STEP 1: Chain A superposition (AF3 predicted fold vs. crystal fold)")
print("=" * 70)
print(f"  RMSD before refinement : {rmsd_before:.3f} A over {n_atoms_before} atoms")
print(f"  RMSD after refinement  : {rmsd_after:.3f} A over {n_atoms_after} atoms")
print(f"  Aligned residues       : {n_res_aligned}")
if rmsd_after < 2.0:
    print(f"  -> PASS sanity check: RMSD {rmsd_after:.3f} A is < 2.0 A; AF3 predicted the fold correctly.")
else:
    print(f"  -> WARNING: RMSD {rmsd_after:.3f} A is >= 2.0 A; fold-level agreement worse than expected.")

# --------------------------------------------------------------------------
# 3. Distance between AF3-predicted phosphate P and crystallographic
#    phosphate P, now that both are in the crystal's reference frame.
# --------------------------------------------------------------------------
crystal_p_sel = find_atom(
    f"crystal and chain {CRYSTAL_LIG_CHAIN} and resi {CRYSTAL_LIG_RESI} and resn {LIG_RESN}",
    P_ATOM_CANDIDATES,
    "crystal PTR phosphate P",
)
af3_p_sel = find_atom(
    f"af3_model and chain {AF3_PEPTIDE_CHAIN} and resi {AF3_LIG_RESI} and resn {LIG_RESN}",
    P_ATOM_CANDIDATES,
    "AF3 predicted PTR phosphate P",
)

p_distance = cmd.get_distance(af3_p_sel, crystal_p_sel, state=1)

print("\n" + "=" * 70)
print("STEP 2: Phosphate (P atom) placement -- crystal PTR vs. AF3-predicted PTR")
print("=" * 70)
print(f"  Distance between phosphate P atoms: {p_distance:.2f} A")
if p_distance < 4.0:
    print("  -> AF3 placed the phosphate group essentially in the same location")
    print("     as the experimentally observed free phosphotyrosine.")
elif p_distance < 8.0:
    print("  -> AF3 placed the phosphate group in the general vicinity of the")
    print("     experimental site, but with a noticeable offset -- inspect visually.")
else:
    print("  -> AF3 placed the phosphate group FAR from the experimentally observed")
    print("     site. It likely did NOT engage the real phosphate-binding pocket.")

# Draw the measured distance in the PyMOL viewer as well.
cmd.distance("ptr_phosphate_distance", af3_p_sel, crystal_p_sel)

# --------------------------------------------------------------------------
# 4. Define the real pocket: crystal protein residues within POCKET_CUTOFF
#    of the crystallographic PTR 301, then test whether the AF3-predicted
#    PTR (post-superposition) contacts those SAME residues.
# --------------------------------------------------------------------------
crystal_lig_sel = f"crystal and chain {CRYSTAL_LIG_CHAIN} and resi {CRYSTAL_LIG_RESI} and resn {LIG_RESN}"
af3_lig_sel = f"af3_model and chain {AF3_PEPTIDE_CHAIN} and resi {AF3_LIG_RESI} and resn {LIG_RESN}"

cmd.select(
    "pocket_residues",
    f"byres (crystal and polymer.protein within {POCKET_CUTOFF} of ({crystal_lig_sel}))",
)
n_pocket_atoms = cmd.count_atoms("pocket_residues")

pocket_list = []
cmd.iterate(
    "pocket_residues and name CA",
    "pocket_list.append((chain, int(resi), resn))",
    space={"pocket_list": pocket_list},
)
pocket_list = sorted(set(pocket_list), key=lambda x: x[1])

print("\n" + "=" * 70)
print(f"STEP 3: Real binding pocket -- crystal residues within {POCKET_CUTOFF} A of PTR 301")
print("=" * 70)
print(f"  {len(pocket_list)} pocket residues identified:")
for chain, resi, resn in pocket_list:
    tag = "  <-- ENGINEERED SUPERBINDER MUTATION" if str(resi) in MUTATED_POSITIONS else ""
    print(f"    {resn} {chain}/{resi}{tag}")

found_mut_positions = {str(r) for _, r, _ in pocket_list}
for mp in MUTATED_POSITIONS:
    status = "IN POCKET (contacts the ligand)" if mp in found_mut_positions else "NOT within cutoff of ligand"
    print(f"  Mutation position {mp}: {status}")

# Which of those SAME pocket residues does the AF3-predicted PTR also contact?
cmd.select(
    "af3_engaged_atoms",
    f"(pocket_residues) within {POCKET_CUTOFF} of ({af3_lig_sel})",
)
cmd.select("af3_engaged_residues", "byres af3_engaged_atoms")

engaged_list = []
cmd.iterate(
    "af3_engaged_residues and name CA",
    "engaged_list.append((chain, int(resi), resn))",
    space={"engaged_list": engaged_list},
)
engaged_list = sorted(set(engaged_list), key=lambda x: x[1])

print("\n" + "=" * 70)
print(f"STEP 4: Does the AF3-predicted PTR also contact those same pocket residues?")
print("=" * 70)
print(f"  {len(engaged_list)} / {len(pocket_list)} real pocket residues are within "
      f"{POCKET_CUTOFF} A of the AF3-predicted PTR:")
for chain, resi, resn in engaged_list:
    tag = "  <-- ENGINEERED SUPERBINDER MUTATION" if str(resi) in MUTATED_POSITIONS else ""
    print(f"    {resn} {chain}/{resi}{tag}")

if pocket_list:
    frac = len(engaged_list) / len(pocket_list)
else:
    frac = 0.0

print(f"\n  Pocket-residue overlap fraction: {frac:.0%}")
if frac >= 0.5 and p_distance < 6.0:
    print("  -> CONCLUSION: AF3's predicted peptide phosphate engages the SAME real")
    print("     superbinder pocket seen crystallographically (good agreement).")
elif frac > 0:
    print("  -> CONCLUSION: AF3's predicted peptide phosphate partially overlaps the")
    print("     real pocket -- inspect the rendered image / session to judge fit.")
else:
    print("  -> CONCLUSION: AF3's predicted peptide phosphate does NOT contact the")
    print("     real superbinder pocket residues -- likely a generic/incorrect pose.")

# --------------------------------------------------------------------------
# 5. Visualization
# --------------------------------------------------------------------------
cmd.hide("everything")
cmd.bg_color("white")
cmd.set("cartoon_transparency", 0.15)

# Crystal: grey cartoon, protein backbone reference
cmd.show("cartoon", "crystal and polymer.protein")
cmd.color("grey80", "crystal and polymer.protein")

# AF3 model: cyan chain A (protein), orange chain B (predicted peptide)
cmd.show("cartoon", f"af3_model and chain {AF3_PROT_CHAIN} and polymer.protein")
cmd.color("cyan", f"af3_model and chain {AF3_PROT_CHAIN} and polymer.protein")
cmd.show("cartoon", f"af3_model and chain {AF3_PEPTIDE_CHAIN} and polymer.protein")
cmd.color("orange", f"af3_model and chain {AF3_PEPTIDE_CHAIN} and polymer.protein")

# Ligands as sticks, distinctly colored
cmd.show("sticks", crystal_lig_sel)
cmd.color("yellow", crystal_lig_sel)
cmd.util.cnc(crystal_lig_sel)

cmd.show("sticks", af3_lig_sel)
cmd.color("magenta", af3_lig_sel)
cmd.util.cnc(af3_lig_sel)

# Real pocket residues as thin sticks/lines for context, and highlight the
# three engineered superbinder positions explicitly.
cmd.show("sticks", "pocket_residues and not name C+N+O")
cmd.color("grey60", "pocket_residues and elem C")
mut_sel = " or ".join(f"(pocket_residues and resi {mp})" for mp in MUTATED_POSITIONS)
cmd.select("mutated_pocket_residues", mut_sel)
cmd.show("sticks", "mutated_pocket_residues and not name C+N+O")
cmd.color("salmon", "mutated_pocket_residues and elem C")

cmd.set("stick_radius", 0.18)
cmd.set("dash_color", "black")

cmd.orient(f"({crystal_lig_sel}) or ({af3_lig_sel})")
cmd.zoom(f"({crystal_lig_sel}) or ({af3_lig_sel}) or pocket_residues", buffer=3)

# --------------------------------------------------------------------------
# 6. Save session + rendered image
# --------------------------------------------------------------------------
os.makedirs(OUTDIR, exist_ok=True)

cmd.save(SESSION_OUT)
print(f"\n[INFO] Session saved: {SESSION_OUT}")

cmd.ray(2400, 2400)
cmd.png(PNG_OUT, dpi=300)
print(f"[INFO] Rendered image saved: {PNG_OUT}")

print("\n" + "=" * 70)
print("DONE.")
print("=" * 70)

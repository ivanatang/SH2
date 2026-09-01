"""
Direct cross-tool comparison of the 8 MD seed structures (4 AF3 medoids +
4 Boltz-2 medoids): superpose all 8 onto a common protein frame, then
compute pairwise peptide backbone RMSD (full peptide and pTyr-anchoring
core Q3-I7) and pTyr side-chain RMSD across all pairs, to see whether the
two tools' representative poses actually agree or occupy different
territory.

Usage: pymol -cq structures/compare_af3_boltz_seeds.py
"""
import json

import numpy as np
from pymol import cmd

BASE = "/Users/ivanatang/Developer/SH2"
SCRATCH = "/private/tmp/claude-501/-Users-ivanatang-Developer-SH2/5b902a39-b797-46c2-89ca-6a1d8551208e/scratchpad"

SEEDS = [
    ("AF3_c1_n85",  f"{BASE}/af3/superB_pYEEI/md_seeds/cluster1_n85_seed-8_sample-2.cif"),
    ("AF3_c2_n53",  f"{BASE}/af3/superB_pYEEI/md_seeds/cluster2_n53_seed-7_sample-2.cif"),
    ("AF3_c3_n7",   f"{BASE}/af3/superB_pYEEI/md_seeds/cluster3_n7_seed-23_sample-4.cif"),
    ("AF3_c4_n5",   f"{BASE}/af3/superB_pYEEI/md_seeds/cluster4_n5_seed-24_sample-2.cif"),
    ("BZ_c1_n133",  f"{BASE}/boltz2/superB_pYEEI/md_seeds/cluster1_n133_model_86.cif"),
    ("BZ_c2_n9",    f"{BASE}/boltz2/superB_pYEEI/md_seeds/cluster2_n9_model_83.cif"),
    ("BZ_c3_n5",    f"{BASE}/boltz2/superB_pYEEI/md_seeds/cluster3_n5_model_111.cif"),
    ("BZ_c4_n3",    f"{BASE}/boltz2/superB_pYEEI/md_seeds/cluster4_n3_model_137.cif"),
]

PTR_SIDECHAIN_ATOMS = "CB+CG+CD1+CD2+CE1+CE2+CZ+OH+P+O1P+O2P+O3P"

cmd.reinitialize()
ref = "AF3_c1_n85"
cmd.load(SEEDS[0][1], ref)

names = []
pep_coords = []
ptr_coords = []
prot_rmsd_to_ref = []

for label, path in SEEDS:
    obj = label
    if obj != ref:
        cmd.load(path, obj)
        aln = cmd.align(f"{obj} and chain A and name N+CA+C+O", f"{ref} and chain A and name N+CA+C+O", cycles=5)
        prot_rmsd_to_ref.append(aln[0])
    else:
        prot_rmsd_to_ref.append(0.0)

    pep_bb_sel = f"{obj} and chain B and name N+CA+C+O and polymer"
    model = cmd.get_model(pep_bb_sel)
    atoms = sorted(model.atom, key=lambda a: (int(a.resi), a.name))
    pep_coords.append([list(a.coord) for a in atoms])

    ptr_sel = f"{obj} and chain B and resn PTR and name {PTR_SIDECHAIN_ATOMS}"
    ptr_model = cmd.get_model(ptr_sel)
    ptr_atoms = sorted(ptr_model.atom, key=lambda a: a.name)
    ptr_coords.append([list(a.coord) for a in ptr_atoms])

    names.append(label)

pep_coords = np.array(pep_coords)   # (8, 44, 3)
ptr_coords = np.array(ptr_coords)   # (8, 12, 3)
n = len(names)

print("=== protein backbone RMSD to reference (AF3_c1_n85), all 8 seeds ===")
for name, r in zip(names, prot_rmsd_to_ref):
    print(f"  {name:12s} {r:.3f} A")

def pairwise(coords):
    m = np.zeros((n, n))
    for i in range(n):
        diff = coords[i][None] - coords
        m[i] = np.sqrt((diff ** 2).sum(axis=2).mean(axis=1))
    return m

mat_full = pairwise(pep_coords)
core_idx = list(range(2 * 4, 7 * 4))  # residues 3-7 (Q3-I7)
mat_core = pairwise(pep_coords[:, core_idx, :])
mat_ptr = pairwise(ptr_coords)

def print_matrix(mat, title):
    print(f"\n=== {title} ===")
    header = "            " + "".join(f"{n:>12s}" for n in names)
    print(header)
    for i, name in enumerate(names):
        row = "".join(f"{mat[i,j]:12.2f}" for j in range(n))
        print(f"{name:12s}{row}")

print_matrix(mat_full, "full peptide backbone RMSD (A)")
print_matrix(mat_core, "core Q3-I7 backbone RMSD (A)")
print_matrix(mat_ptr, "pTyr side-chain RMSD (A)")

af3_idx = [i for i, n_ in enumerate(names) if n_.startswith("AF3")]
bz_idx = [i for i, n_ in enumerate(names) if n_.startswith("BZ")]

def block_stats(mat, rows, cols, symmetric=False):
    vals = []
    for i in rows:
        for j in cols:
            if symmetric and i >= j:
                continue
            if not symmetric and i == j:
                continue
            vals.append(mat[i, j])
    vals = np.array(vals)
    return vals.mean(), vals.min(), vals.max()

print("\n=== summary: within-AF3 vs within-Boltz2 vs cross-tool ===")
for label, mat in [("full peptide", mat_full), ("core Q3-I7", mat_core), ("pTyr sidechain", mat_ptr)]:
    wa = block_stats(mat, af3_idx, af3_idx, symmetric=True)
    wb = block_stats(mat, bz_idx, bz_idx, symmetric=True)
    cr = block_stats(mat, af3_idx, bz_idx, symmetric=False)
    print(f"{label:16s}  within-AF3: mean={wa[0]:.2f} [{wa[1]:.2f},{wa[2]:.2f}]   "
          f"within-Boltz2: mean={wb[0]:.2f} [{wb[1]:.2f},{wb[2]:.2f}]   "
          f"cross-tool: mean={cr[0]:.2f} [{cr[1]:.2f},{cr[2]:.2f}]")

np.savez(f"{SCRATCH}/seed_comparison.npz", mat_full=mat_full, mat_core=mat_core, mat_ptr=mat_ptr)
with open(f"{SCRATCH}/seed_comparison_names.json", "w") as fh:
    json.dump(names, fh)

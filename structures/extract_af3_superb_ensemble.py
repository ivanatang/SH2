"""
Superpose all 150 AF3 superbinder-SH2 + pYEEI models onto a common protein
frame (chain A backbone), then dump the peptide (chain B) backbone + PTR
side-chain coordinates plus protein self-consistency RMSD to JSON, for
downstream clustering in a normal Python env (pymol's bundled python has no
scipy/sklearn).

Usage: pymol -cq structures/extract_af3_superb_ensemble.py
"""
import glob
import json
import os

import numpy as np
from pymol import cmd

BASE = "/Users/ivanatang/Developer/SH2"
OUT_DIR = f"{BASE}/af3/superB_pYEEI/output/superbinder_sh2_epqyeeipiyl"
SCRATCH = "/private/tmp/claude-501/-Users-ivanatang-Developer-SH2/5b902a39-b797-46c2-89ca-6a1d8551208e/scratchpad"
OUT_JSON = f"{SCRATCH}/af3_superb_ensemble.json"

model_dirs = sorted(glob.glob(f"{OUT_DIR}/seed-*_sample-*"))
assert len(model_dirs) == 150, f"expected 150, found {len(model_dirs)}"

PTR_SIDECHAIN_ATOMS = "CB+CG+CD1+CD2+CE1+CE2+CZ+OH+P+O1P+O2P+O3P"

cmd.reinitialize()
ref_name = "ref"
cmd.load(f"{model_dirs[0]}/model.cif", ref_name)
cmd.remove(f"{ref_name} and not chain A")

results = []
for d in model_dirs:
    name = os.path.basename(d)
    obj = "m_" + name.replace("-", "_")
    cmd.load(f"{d}/model.cif", obj)

    aln = cmd.align(
        f"{obj} and chain A and name N+CA+C+O",
        f"{ref_name} and name N+CA+C+O",
        cycles=5,
    )
    protein_rmsd_after = aln[0]

    pep_bb_sel = f"{obj} and chain B and name N+CA+C+O and polymer"
    model = cmd.get_model(pep_bb_sel)
    atoms = sorted(model.atom, key=lambda a: (int(a.resi), a.name))
    pep_coords = [list(a.coord) for a in atoms]

    ptr_sel = f"{obj} and chain B and resn PTR and name {PTR_SIDECHAIN_ATOMS}"
    ptr_model = cmd.get_model(ptr_sel)
    ptr_atoms = sorted(ptr_model.atom, key=lambda a: a.name)
    ptr_coords = [list(a.coord) for a in ptr_atoms]

    results.append(
        dict(
            name=name,
            protein_rmsd_to_ref=protein_rmsd_after,
            pep_bb_coords=pep_coords,
            ptr_sc_coords=ptr_coords,
        )
    )
    cmd.delete(obj)

prot_rmsds = np.array([r["protein_rmsd_to_ref"] for r in results])
print(f"protein backbone RMSD to reference (seed-1_sample-0), n={len(results)}:")
print(f"  mean={prot_rmsds.mean():.3f} std={prot_rmsds.std():.3f} max={prot_rmsds.max():.3f}")

with open(OUT_JSON, "w") as fh:
    json.dump(results, fh)
print(f"saved -> {OUT_JSON}")

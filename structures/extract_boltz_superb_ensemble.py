"""
Same pipeline as extract_af3_superb_ensemble.py, applied to the 150 Boltz-2
superbinder-SH2 + pYEEI samples (boltz2/superB_pYEEI/output/..., seed=1,
diffusion_samples=150). Superposes all 150 onto a common protein (chain A)
frame and dumps peptide (chain B) backbone + PTR side-chain coordinates to
JSON for downstream clustering.

Usage: pymol -cq structures/extract_boltz_superb_ensemble.py
"""
import json
import os

import numpy as np
from pymol import cmd

BASE = "/Users/ivanatang/Developer/SH2"
PRED_DIR = (
    f"{BASE}/boltz2/superB_pYEEI/output/boltz_results_superbinder_sh2_pYEEI/"
    f"predictions/superbinder_sh2_pYEEI"
)
OUT_JSON = f"{BASE}/structures/boltz_superb_ensemble.json"

N_MODELS = 150
PTR_SIDECHAIN_ATOMS = "CB+CG+CD1+CD2+CE1+CE2+CZ+OH+P+O1P+O2P+O3P"

cmd.reinitialize()
ref_name = "ref"
cmd.load(f"{PRED_DIR}/superbinder_sh2_pYEEI_model_0.cif", ref_name)
cmd.remove(f"{ref_name} and not chain A")

results = []
for i in range(N_MODELS):
    name = f"model_{i}"
    obj = f"m_{i}"
    cmd.load(f"{PRED_DIR}/superbinder_sh2_pYEEI_model_{i}.cif", obj)

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
print(f"protein backbone RMSD to reference (model_0), n={len(results)}:")
print(f"  mean={prot_rmsds.mean():.3f} std={prot_rmsds.std():.3f} max={prot_rmsds.max():.3f}")

with open(OUT_JSON, "w") as fh:
    json.dump(results, fh)
print(f"saved -> {OUT_JSON}")

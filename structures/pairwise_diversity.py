"""
Pairwise peptide-pose diversity within the AF3 ensemble (5 models) and within
the Boltz-2 ensemble (25 samples): after superposing every model onto a common
protein frame (AF3 model_0's SH2 chain), how much does the peptide's pose vary
model-to-model within each tool's own output set? This tells us whether either
tool is actually proposing distinct alternative binding poses (useful as
diverse MD seeds) or just re-generating the same pose with noise.

Usage: pymol -cq structures/pairwise_diversity.py
"""
import itertools
import json

import numpy as np
from pymol import cmd

BASE = "/Users/ivanatang/Developer/SH2"
AF3_DIR = f"{BASE}/af3/fold_wt_sh2_pyeei"
BOLTZ_DIR = (
    f"{BASE}/boltz2/WT_SH2_pYEEI/boltz_cSrc_validation/"
    f"boltz_results_cSrc_SH2_EPQpYEEIPIYL/predictions/cSrc_SH2_EPQpYEEIPIYL"
)

cmd.reinitialize()
cmd.load(f"{AF3_DIR}/fold_wt_sh2_pyeei_model_0.cif", "ref_frame")
cmd.remove("ref_frame and not chain A")

entries = []  # (label, obj_name)

for i in range(5):
    name = f"af3_{i}"
    cmd.load(f"{AF3_DIR}/fold_wt_sh2_pyeei_model_{i}.cif", name)
    cmd.align(f"{name} and chain A and name N+CA+C+O", "ref_frame and name N+CA+C+O", cycles=5)
    entries.append((f"AF3_m{i}", name))

for i in range(25):
    name = f"bz_{i}"
    cmd.load(f"{BOLTZ_DIR}/cSrc_SH2_EPQpYEEIPIYL_model_{i}.cif", name)
    cmd.align(f"{name} and chain A and name N+CA+C+O", "ref_frame and name N+CA+C+O", cycles=5)
    entries.append((f"BZ_m{i}", name))

# full-length peptide backbone (all 11 residues -- predictions are fully resolved)
def get_coords(obj):
    sel = f"{obj} and chain B and name N+CA+C+O and polymer"
    model = cmd.get_model(sel)
    atoms = sorted(model.atom, key=lambda a: (int(a.resi), a.name))
    return np.array([a.coord for a in atoms])

coords = {label: get_coords(obj) for label, obj in entries}
n = len(entries)
labels = [l for l, _ in entries]
mat = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i == j:
            continue
        a, b = coords[labels[i]], coords[labels[j]]
        m = min(len(a), len(b))
        mat[i, j] = np.sqrt(((a[:m] - b[:m]) ** 2).sum(axis=1).mean())

af3_idx = [k for k, l in enumerate(labels) if l.startswith("AF3")]
bz_idx = [k for k, l in enumerate(labels) if l.startswith("BZ")]


def offdiag_stats(idx):
    vals = [mat[i, j] for i in idx for j in idx if i != j]
    return np.mean(vals), np.std(vals), np.min(vals), np.max(vals)


def cross_stats(idx1, idx2):
    vals = [mat[i, j] for i in idx1 for j in idx2]
    return np.mean(vals), np.min(vals), np.max(vals)


print("=== within-AF3 pairwise peptide backbone RMSD (5 models, common protein frame) ===")
m, s, lo, hi = offdiag_stats(af3_idx)
print(f"mean={m:.3f} std={s:.3f} min={lo:.3f} max={hi:.3f}")

print("\n=== within-Boltz2 pairwise peptide backbone RMSD (25 models, common protein frame) ===")
m, s, lo, hi = offdiag_stats(bz_idx)
print(f"mean={m:.3f} std={s:.3f} min={lo:.3f} max={hi:.3f}")

print("\n=== AF3 vs Boltz2 cross-set pairwise peptide backbone RMSD ===")
m, lo, hi = cross_stats(af3_idx, bz_idx)
print(f"mean={m:.3f} min={lo:.3f} max={hi:.3f}")

with open(f"{BASE}/structures/pairwise_rmsd_matrix.json", "w") as fh:
    json.dump(dict(labels=labels, matrix=mat.tolist()), fh)
print(f"\nSaved pairwise matrix -> {BASE}/structures/pairwise_rmsd_matrix.json")

"""
Cluster the 150 AF3 superbinder-SH2 + pYEEI models (30 seeds x 5 samples,
af3/superB_pYEEI/output/superbinder_sh2_epqyeeipiyl/) by peptide binding pose,
and pick a medoid representative per cluster to use as an MD seed structure.

Two-step pipeline (pymol has no scipy/sklearn, so coordinate extraction and
clustering are separate scripts):

  1. structures/extract_af3_superb_ensemble.py  (run with `pymol -cq`)
     Superposes all 150 models onto a common protein (chain A) frame and
     dumps peptide (chain B) backbone + PTR side-chain coordinates to JSON.

  2. This script (run with a normal python env that has numpy/scipy/sklearn)
     Builds the pairwise peptide-backbone RMSD matrix, does average-linkage
     hierarchical clustering, and reports medoids for several k.

Key finding from the WT-vs-1SPS validation carried over here: cluster on the
peptide pose, not the protein (protein backbone was self-consistent to
~0.13 A across all 150 -- essentially zero protein-side diversity to
cluster on). Also worth knowing before picking k: the pTyr-anchoring
specificity core (residues Q3-I7) is essentially invariant across all 150
samples (pairwise RMSD <= 0.86 A, none >1.0 A) -- clustering signal comes
entirely from the flexible N-/C-terminal tails (E1-P2 dangling end
especially), which is also the region 1SPS itself left unresolved in the
crystal. So more clusters = more coverage of terminus orientation, not
evidence of alternate core binding modes.

Usage:
    python3 structures/cluster_af3_superb_ensemble.py [--k K]
"""
import argparse
import json

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score

BASE = "/Users/ivanatang/Developer/SH2"
ENSEMBLE_JSON = f"{BASE}/structures/af3_superb_ensemble.json"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=None, help="number of clusters; default sweeps 2-15 and reports silhouette")
    args = parser.parse_args()

    with open(ENSEMBLE_JSON) as fh:
        data = json.load(fh)
    names = [d["name"] for d in data]
    n = len(names)
    pep_coords = np.array([d["pep_bb_coords"] for d in data])  # (n, 44, 3): 11 residues x N,CA,C,O

    mat = np.zeros((n, n))
    for i in range(n):
        diff = pep_coords[i][None, :, :] - pep_coords
        mat[i] = np.sqrt((diff ** 2).sum(axis=2).mean(axis=1))

    condensed = squareform(mat, checks=False)
    Z = linkage(condensed, method="average")

    k_values = [args.k] if args.k else range(2, 16)
    best_k, best_score = None, -2
    for k in k_values:
        labels = fcluster(Z, t=k, criterion="maxclust")
        if len(set(labels)) < 2:
            continue
        score = silhouette_score(mat, labels, metric="precomputed")
        sizes = sorted(np.bincount(labels)[1:], reverse=True)
        print(f"k={k:2d}  silhouette={score:.3f}  sizes={sizes}")
        if score > best_score:
            best_score, best_k = score, k

    k = args.k or best_k
    labels = fcluster(Z, t=k, criterion="maxclust")
    clusters = {}
    for lbl, name in zip(labels, names):
        clusters.setdefault(int(lbl), []).append(name)

    medoids = {}
    for lbl, members in clusters.items():
        idx = [names.index(m) for m in members]
        sub = mat[np.ix_(idx, idx)]
        medoid_local = idx[np.argmin(sub.sum(axis=1))]
        medoids[lbl] = names[medoid_local]

    print(f"\n=== final clustering, k={k} ===")
    for lbl in sorted(clusters, key=lambda l: -len(clusters[l])):
        members = clusters[lbl]
        idx = [names.index(m) for m in members]
        sub = mat[np.ix_(idx, idx)]
        intra = sub[np.triu_indices(len(idx), k=1)].mean() if len(idx) > 1 else 0.0
        print(f"  cluster {lbl}: n={len(members):3d} ({100*len(members)/n:.1f}%)  intra-mean={intra:.2f}A  medoid={medoids[lbl]}")

    out = dict(
        k=k,
        clusters={str(l): {"members": m, "medoid": medoids[l]} for l, m in clusters.items()},
    )
    out_path = f"{BASE}/structures/af3_superb_clusters_k{k}.json"
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nsaved -> {out_path}")


if __name__ == "__main__":
    main()

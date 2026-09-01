"""
Cluster the 150 Boltz-2 superbinder-SH2 + pYEEI models (seed=1,
diffusion_samples=150, boltz2/superB_pYEEI/output/...) by peptide binding
pose, and pick the medoid of each cluster as an MD seed structure.

Same two-step pipeline as the AF3 version
(structures/cluster_af3_superb_ensemble.py):
  1. structures/extract_boltz_superb_ensemble.py (run with `pymol -cq`)
     Superposes all 150 models onto a common protein frame, dumps peptide
     backbone + PTR side-chain coordinates to JSON.
  2. This script (normal python env with numpy/scipy/sklearn)
     Pairwise peptide-backbone RMSD, average-linkage clustering, medoids.

Usage:
    python3 structures/cluster_boltz_superb_ensemble.py [--k K]
"""
import argparse
import json

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score

BASE = "/Users/ivanatang/Developer/SH2"
ENSEMBLE_JSON = f"{BASE}/structures/boltz_superb_ensemble.json"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=None)
    args = parser.parse_args()

    with open(ENSEMBLE_JSON) as fh:
        data = json.load(fh)
    names = [d["name"] for d in data]
    n = len(names)
    pep_coords = np.array([d["pep_bb_coords"] for d in data])

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

    out = dict(k=k, clusters={str(l): {"members": m, "medoid": medoids[l]} for l, m in clusters.items()})
    out_path = f"{BASE}/structures/boltz_superb_clusters_k{k}.json"
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nsaved -> {out_path}")


if __name__ == "__main__":
    main()

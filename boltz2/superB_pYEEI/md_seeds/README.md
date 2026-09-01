# Boltz-2 superbinder-SH2 + pYEEI: MD seed structures

4 representative binding conformations, selected by clustering all 150
Boltz-2 predictions (`../superbinder_sh2_pYEEI.yaml`, seed=1,
diffusion_samples=150, `run_boltz.sh`) and picking the medoid of each
cluster. Same methodology as `../../../af3/superB_pYEEI/md_seeds/README.md`
-- scripts: `../../../structures/extract_boltz_superb_ensemble.py`
(coordinate extraction/alignment) and
`../../../structures/cluster_boltz_superb_ensemble.py` (clustering).

## QC on the 150 predictions

ptm 0.95-0.98, iptm 0.90-0.96 (mean 0.93), no outliers. Note: this run's
confidence JSON schema doesn't include `has_clash`/`fraction_disordered`
(different Boltz-2 version than the earlier WT validation run), so those
specific checks weren't available here -- ptm/iptm gave no reason to
flag or discard any sample. Protein backbone is self-consistent across
the ensemble (mean RMSD to a reference sample: 0.13 +/- 0.03 A),
matching the AF3 result almost exactly -- clustering was done on the
peptide pose only, as before.

## Same core finding as AF3, opposite tail

The pTyr-anchoring specificity core (residues Q3-I7) is essentially
invariant across all 150 samples -- pairwise RMSD never exceeds 0.67 A,
0% of pairs exceed 1.0 A (even tighter than AF3's 0.86 A max). All AF3
and Boltz-2 samples agree on the core bound pose.

Where the two tools disagree is *which* flexible terminus carries the
diversity. AF3's ensemble varied mainly at the N-terminal Glu1-Pro2
dangling end (up to ~6 A swings). Boltz-2's ensemble instead varies
mainly at the **C-terminal Ile9-Tyr10-Leu11** tail (per-residue CA
std: Leu11 2.41 A, Tyr10 1.30 A -- vs. Glu1 only 0.79 A here). Both
termini are regions 1SPS itself couldn't fully resolve, so both are
plausible real flexibility; the two tools are independently sampling
different parts of it rather than disagreeing about the core.

## Clustering

Average-linkage hierarchical clustering on peptide backbone RMSD (all 11
residues, common protein frame), same as the AF3 pipeline. Silhouette
again peaks at k=2 (0.45), but the split is far more lopsided than AF3's
(141/9 vs. AF3's balanced 92/58). k=4 was used to match the AF3 pipeline's
seed count for comparability rather than re-deriving k from scratch here;
full silhouette-vs-k sweep and cluster membership at k=2/3/4/5 (plus the
statistically-best k):
`../../../structures/boltz_superb_clusters_multi_k.json`.

| file | cluster | size (of 150) | intra-cluster mean RMSD |
|---|---|---|---|
| `cluster1_n133_model_86.cif` | 1 (dominant) | 133 (88.7%) | 1.15 A |
| `cluster2_n9_model_83.cif` | 2 (minor) | 9 (6.0%) | 1.12 A |
| `cluster3_n5_model_111.cif` | 3 (minor) | 5 (3.3%) | 1.09 A |
| `cluster4_n3_model_137.cif` | 4 (minor) | 3 (2.0%) | 0.85 A |

Unlike AF3's split (two comparably-sized populations), Boltz-2's k=4
clustering is dominated by one large cluster (88.7%) with three small
offshoots -- worth keeping in mind when weighting these 4 seeds for MD
(cluster 1's terminus orientation is far more representative of the
Boltz-2 ensemble than clusters 2-4 are).

Each file is the medoid (minimum summed RMSD to every other member of its
cluster) -- an actual predicted structure, not an averaged/synthetic one.

## Next step

Same as the AF3 seeds: feed each of these 4 structures through the
`md_prep/` system-prep pipeline to set up independent MD systems. With
both AF3 and Boltz-2 seed sets now in hand (8 structures total across the
two tools), worth deciding whether to run MD from all 8 or select a
subset -- e.g. both tools' dominant-cluster medoid plus 1-2 minor-cluster
representatives from each, to cover both the N-terminal (AF3) and
C-terminal (Boltz-2) flexibility observed.

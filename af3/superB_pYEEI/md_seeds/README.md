# AF3 superbinder-SH2 + pYEEI: MD seed structures

4 representative binding conformations, selected by clustering all 150 AF3
predictions (`../superbinder_sh2_pYEEI.json`, 30 seeds x 5 samples,
`run_af3_sh2.sh`) and picking the medoid of each cluster. Full methodology
and scripts: `../../../structures/extract_af3_superb_ensemble.py` (coordinate
extraction/alignment) and `../../../structures/cluster_af3_superb_ensemble.py`
(clustering).

## QC on the 150 predictions

All 150 passed quality screening before clustering: `has_clash = 0` for
every sample, ptm 0.88-0.90, iptm 0.86-0.89 (tight distributions, no
outliers to discard). Protein backbone is self-consistent across the whole
ensemble (mean RMSD to a reference sample: 0.13 +/- 0.02 A) -- essentially
no protein-side conformational diversity, so clustering was done on the
peptide binding pose only, after superposing every model onto a common
protein frame.

## Key finding: clustering signal is terminus flexibility, not pocket ambiguity

The pTyr-anchoring specificity core (peptide residues Q3-I7: the pTyr
pocket + the immediate +1/+2/+3 contacts that define SH2 binding
specificity) is essentially identical across all 150 samples -- pairwise
RMSD in that region never exceeds 0.86 A. All the clustering signal (mean
pairwise peptide backbone RMSD 1.94 A across the full ensemble) comes from
the flexible N-terminal Glu1-Pro2 dangling end (up to ~6 A swings) and, to
a lesser extent, the C-terminal Ile9-Tyr10-Leu11 tail. This is the same
region 1SPS (the WT crystal reference used in the AF3-vs-Boltz2 validation)
itself left unresolved/disordered. So: AF3 unanimously agrees on the core
bound pose: the 4 clusters below differ only in how the flexible termini
are oriented, which is legitimate structural diversity worth seeding MD
from (rather than letting every replica start with the tail in one
arbitrary orientation), not evidence of a genuinely ambiguous binding mode.

## Clustering

Average-linkage hierarchical clustering on peptide backbone RMSD (all 11
residues, common protein frame). Silhouette score peaks at k=2 (0.50); k=4
was chosen to also capture two smaller, less dominant terminus populations
before cluster quality degrades into noise (silhouette keeps falling past
k=4). Full silhouette-vs-k sweep and cluster membership at k=2/3/4/5:
`../../../structures/af3_superb_clusters_multi_k.json`.

| file | cluster | size (of 150) | intra-cluster mean RMSD |
|---|---|---|---|
| `cluster1_n85_seed-8_sample-2.cif` | 1 (largest) | 85 (56.7%) | 1.15 A |
| `cluster2_n53_seed-7_sample-2.cif` | 2 | 53 (35.3%) | 1.31 A |
| `cluster3_n7_seed-23_sample-4.cif` | 3 (minor) | 7 (4.7%) | 1.18 A |
| `cluster4_n5_seed-24_sample-2.cif` | 4 (minor) | 5 (3.3%) | 1.38 A |

Each file is the medoid (the ensemble member with minimum summed RMSD to
every other member of its cluster) -- an actual predicted structure, not an
averaged/synthetic one, so it's directly usable as an MD starting
coordinate set.

## Next step

Feed each of these 4 structures through the existing single-structure prep
pipeline in `md_prep/` (parameterization, solvation, EM/equil/production)
to set up 4 independent MD systems. The Boltz-2 150-sample run
(`boltz2/superB_pYEEI/`) is still in progress; once it finishes, the same
two-script pipeline can be applied there and the results compared/merged
with this AF3 clustering.

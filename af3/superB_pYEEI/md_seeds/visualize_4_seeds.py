"""
Visualize the 4 cluster-medoid MD seed structures overlaid on a common
protein frame, so peptide pose differences between clusters are clearly
visible (and the shared, invariant pTyr pocket engagement is too).

Usage: pymol -cq af3/superB_pYEEI/md_seeds/visualize_4_seeds.py
"""
from pymol import cmd

BASE = "/Users/ivanatang/Developer/SH2"
SEED_DIR = f"{BASE}/af3/superB_pYEEI/md_seeds"
OUT_DIR = SEED_DIR

# (object name, file, color, label)
CLUSTERS = [
    ("c1", f"{SEED_DIR}/cluster1_n85_seed-8_sample-2.cif", "marine", "cluster 1 (n=85, 56.7%)"),
    ("c2", f"{SEED_DIR}/cluster2_n53_seed-7_sample-2.cif", "orange", "cluster 2 (n=53, 35.3%)"),
    ("c3", f"{SEED_DIR}/cluster3_n7_seed-23_sample-4.cif", "forest", "cluster 3 (n=7, 4.7%)"),
    ("c4", f"{SEED_DIR}/cluster4_n5_seed-24_sample-2.cif", "purple", "cluster 4 (n=5, 3.3%)"),
]

cmd.reinitialize()
cmd.bg_color("white")
cmd.set("ray_opaque_background", 0)
cmd.set("cartoon_transparency", 0.0)

for name, path, color, label in CLUSTERS:
    cmd.load(path, name)
    cmd.dss(name)  # assign secondary structure so the protein cartoon shows helices/sheets, not a bare tube

# superpose all 4 onto cluster 1 (the majority pose) protein backbone
ref = "c1"
for name, *_ in CLUSTERS[1:]:
    cmd.align(f"{name} and chain A and name N+CA+C+O", f"{ref} and chain A and name N+CA+C+O", cycles=5)

# protein: show once, from the reference, as a pale transparent cartoon --
# all 4 protein chains are self-consistent to ~0.13 A so overlaying 4 would
# just clutter the image without adding information
cmd.hide("everything")
cmd.show("cartoon", f"{ref} and chain A and polymer.protein")
cmd.color("gray80", f"{ref} and chain A and polymer.protein")
cmd.set("cartoon_transparency", 0.55, f"{ref} and chain A")

# peptides: sticks only (cleaner than cartoon for an 11-mer), one color per
# cluster, PTR phosphate emphasized as spheres since it's the shared anchor
for name, path, color, label in CLUSTERS:
    sel = f"{name} and chain B"
    cmd.show("sticks", f"{sel} and not hydro")
    cmd.color(color, f"{sel} and elem C")
    cmd.util.cnc(sel)
    cmd.show("spheres", f"{sel} and resn PTR and name P+O1P+O2P+O3P")
    cmd.set("sphere_scale", 0.35, f"{sel} and resn PTR and name P+O1P+O2P+O3P")

cmd.orient(f"{ref} and chain A or c1 and chain B")
cmd.zoom("chain A or chain B", buffer=4)
cmd.set("ray_trace_mode", 1)
cmd.ray(2000, 2000)
cmd.png(f"{OUT_DIR}/4clusters_overlay_full.png", dpi=300)

# tight close-up on just the peptide, protein hidden, to show terminus divergence clearly
cmd.hide("everything", "chain A")
cmd.orient("chain B")
cmd.zoom("chain B", buffer=3)
cmd.ray(2000, 2000)
cmd.png(f"{OUT_DIR}/4clusters_overlay_peptide_closeup.png", dpi=300)

# restore protein for the saved session
cmd.show("cartoon", f"{ref} and chain A and polymer.protein")
cmd.save(f"{OUT_DIR}/4clusters_overlay.pse")

print(f"saved -> {OUT_DIR}/4clusters_overlay_full.png")
print(f"saved -> {OUT_DIR}/4clusters_overlay_peptide_closeup.png")
print(f"saved -> {OUT_DIR}/4clusters_overlay.pse")

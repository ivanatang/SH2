"""
Visualize the 4 Boltz-2 cluster-medoid MD seed structures overlaid on a
common protein frame. Mirrors af3/superB_pYEEI/md_seeds/visualize_4_seeds.py.

Usage: pymol -cq boltz2/superB_pYEEI/md_seeds/visualize_4_seeds.py
"""
from pymol import cmd

BASE = "/Users/ivanatang/Developer/SH2"
SEED_DIR = f"{BASE}/boltz2/superB_pYEEI/md_seeds"
OUT_DIR = SEED_DIR

CLUSTERS = [
    ("c1", f"{SEED_DIR}/cluster1_n133_model_86.cif", "marine", "cluster 1 (n=133, 88.7%)"),
    ("c2", f"{SEED_DIR}/cluster2_n9_model_83.cif", "orange", "cluster 2 (n=9, 6.0%)"),
    ("c3", f"{SEED_DIR}/cluster3_n5_model_111.cif", "forest", "cluster 3 (n=5, 3.3%)"),
    ("c4", f"{SEED_DIR}/cluster4_n3_model_137.cif", "purple", "cluster 4 (n=3, 2.0%)"),
]

cmd.reinitialize()
cmd.bg_color("white")
cmd.set("ray_opaque_background", 0)

for name, path, color, label in CLUSTERS:
    cmd.load(path, name)
    cmd.dss(name)

ref = "c1"
for name, *_ in CLUSTERS[1:]:
    cmd.align(f"{name} and chain A and name N+CA+C+O", f"{ref} and chain A and name N+CA+C+O", cycles=5)

cmd.hide("everything")
cmd.show("cartoon", f"{ref} and chain A and polymer.protein")
cmd.color("gray80", f"{ref} and chain A and polymer.protein")
cmd.set("cartoon_transparency", 0.55, f"{ref} and chain A")

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

cmd.hide("everything", "chain A")
cmd.orient("chain B")
cmd.zoom("chain B", buffer=3)
cmd.ray(2000, 2000)
cmd.png(f"{OUT_DIR}/4clusters_overlay_peptide_closeup.png", dpi=300)

cmd.show("cartoon", f"{ref} and chain A and polymer.protein")
cmd.save(f"{OUT_DIR}/4clusters_overlay.pse")

print(f"saved -> {OUT_DIR}/4clusters_overlay_full.png")
print(f"saved -> {OUT_DIR}/4clusters_overlay_peptide_closeup.png")
print(f"saved -> {OUT_DIR}/4clusters_overlay.pse")

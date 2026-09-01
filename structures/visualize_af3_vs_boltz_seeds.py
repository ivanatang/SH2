"""
Combined overlay of all 8 MD seed structures (4 AF3 + 4 Boltz-2), colored by
tool rather than by cluster, to visualize the cross-tool comparison:
core backbone register agreement, terminus-flexibility differences
(AF3=N-term, Boltz2=C-term), and pTyr side-chain rotamer spread.

Usage: pymol -cq structures/visualize_af3_vs_boltz_seeds.py
"""
from pymol import cmd

BASE = "/Users/ivanatang/Developer/SH2"
OUT_DIR = f"{BASE}/structures"

SEEDS = [
    ("AF3_c1", f"{BASE}/af3/superB_pYEEI/md_seeds/cluster1_n85_seed-8_sample-2.cif", "skyblue"),
    ("AF3_c2", f"{BASE}/af3/superB_pYEEI/md_seeds/cluster2_n53_seed-7_sample-2.cif", "skyblue"),
    ("AF3_c3", f"{BASE}/af3/superB_pYEEI/md_seeds/cluster3_n7_seed-23_sample-4.cif", "skyblue"),
    ("AF3_c4", f"{BASE}/af3/superB_pYEEI/md_seeds/cluster4_n5_seed-24_sample-2.cif", "skyblue"),
    ("BZ_c1", f"{BASE}/boltz2/superB_pYEEI/md_seeds/cluster1_n133_model_86.cif", "orange"),
    ("BZ_c2", f"{BASE}/boltz2/superB_pYEEI/md_seeds/cluster2_n9_model_83.cif", "orange"),
    ("BZ_c3", f"{BASE}/boltz2/superB_pYEEI/md_seeds/cluster3_n5_model_111.cif", "orange"),
    ("BZ_c4", f"{BASE}/boltz2/superB_pYEEI/md_seeds/cluster4_n3_model_137.cif", "orange"),
]

cmd.reinitialize()
cmd.bg_color("white")
cmd.set("ray_opaque_background", 0)

ref = "AF3_c1"
for name, path, color in SEEDS:
    cmd.load(path, name)
    cmd.dss(name)
    if name != ref:
        cmd.align(f"{name} and chain A and name N+CA+C+O", f"{ref} and chain A and name N+CA+C+O", cycles=5)

cmd.hide("everything")
cmd.show("cartoon", f"{ref} and chain A and polymer.protein")
cmd.color("gray85", f"{ref} and chain A and polymer.protein")
cmd.set("cartoon_transparency", 0.6, f"{ref} and chain A")

for name, path, color in SEEDS:
    sel = f"{name} and chain B"
    cmd.show("sticks", f"{sel} and not hydro")
    cmd.color(color, f"{sel} and elem C")
    cmd.util.cnc(sel)
    cmd.show("spheres", f"{sel} and resn PTR and name P+O1P+O2P+O3P")
    cmd.set("sphere_scale", 0.35, f"{sel} and resn PTR and name P+O1P+O2P+O3P")
    cmd.set("stick_radius", 0.22, sel)

cmd.orient(f"{ref} and chain A or {ref} and chain B")
cmd.zoom("chain A or chain B", buffer=4)
cmd.set("ray_trace_mode", 1)
cmd.ray(2000, 2000)
cmd.png(f"{OUT_DIR}/af3_vs_boltz_seeds_full.png", dpi=300)

cmd.hide("everything", "chain A")
cmd.orient("chain B")
cmd.zoom("chain B", buffer=3)
cmd.ray(2000, 2000)
cmd.png(f"{OUT_DIR}/af3_vs_boltz_seeds_peptide_closeup.png", dpi=300)

# extra-tight close-up on just the pTyr ring + phosphate, to show rotamer spread
cmd.orient("chain B and resn PTR")
cmd.zoom("chain B and resn PTR", buffer=2.5)
cmd.ray(2000, 2000)
cmd.png(f"{OUT_DIR}/af3_vs_boltz_seeds_ptyr_closeup.png", dpi=300)

cmd.show("cartoon", f"{ref} and chain A and polymer.protein")
cmd.save(f"{OUT_DIR}/af3_vs_boltz_seeds.pse")
print("done")

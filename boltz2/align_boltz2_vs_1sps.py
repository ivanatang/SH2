"""
align_boltz2_vs_1sps.py

Structural validation: how well do the 25 Boltz-2 diffusion samples of WT
c-Src SH2 + phosphopeptide pYEEI (EPQ-pY-EEIPIYL) reproduce the experimental
crystal structure 1SPS (Waksman & Kuriyan 1993, 2.7 A)?

Mirrors af3/align_af3_vs_1sps_5models.py exactly (same reference prep, same
metrics, same peptide numbering convention) so the two model families are
directly comparable.

Inputs
------
- Boltz-2 models (mmCIF), chain A = SH2 protein (104 res), chain B = peptide
  (11 res, PTR at resi 4), auth_asym_id/auth_seq_id verified identical
  convention to the AF3 files:
    boltz2/WT_SH2_pYEEI/boltz_cSrc_validation/boltz_results_cSrc_SH2_EPQpYEEIPIYL/
    predictions/cSrc_SH2_EPQpYEEIPIYL/cSrc_SH2_EPQpYEEIPIYL_model_{0..24}.cif
  (25 diffusion samples from one seed; confidence_score/ptm/iptm read from the
  sibling confidence_..._model_N.json files.)
- Crystal reference: structures/1SPS.pdb (chain A protein / chain D peptide,
  same prep as the AF3 script: peptide chain D renumbered +4 onto the 1-11
  frame, only resolved peptide positions 2-8 used for peptide RMSD).

Usage
-----
    pymol -cq boltz2/align_boltz2_vs_1sps.py

Outputs
-------
    structures/boltz2_vs_1sps_model0_overlay.png / .pse   (top-confidence model)
    structures/boltz2_vs_1sps_results.json                (all 25 models)
"""

import glob
import json
import os

import numpy as np
from pymol import cmd

BASE = "/Users/ivanatang/Developer/SH2"
REF_PDB = f"{BASE}/structures/1SPS.pdb"
BOLTZ_DIR = (
    f"{BASE}/boltz2/WT_SH2_pYEEI/boltz_cSrc_validation/"
    f"boltz_results_cSrc_SH2_EPQpYEEIPIYL/predictions/cSrc_SH2_EPQpYEEIPIYL"
)
OUT_PNG = f"{BASE}/structures/boltz2_vs_1sps_model0_overlay.png"
OUT_PSE = f"{BASE}/structures/boltz2_vs_1sps_model0_overlay.pse"
OUT_JSON = f"{BASE}/structures/boltz2_vs_1sps_results.json"

N_MODELS = 25

PEP_RESI = list(range(2, 9))  # 2..8, same resolved/comparable region as AF3 script
PTR_RESI = 4
PTR_SIDECHAIN_ATOMS = "CB+CG+CD1+CD2+CE1+CE2+CZ+OH+P+O1P+O2P+O3P"

# read Boltz-2's own confidence metrics so we can rank vs. crystallographic accuracy
conf_by_model = {}
for i in range(N_MODELS):
    cpath = f"{BOLTZ_DIR}/confidence_cSrc_SH2_EPQpYEEIPIYL_model_{i}.json"
    with open(cpath) as fh:
        c = json.load(fh)
    conf_by_model[i] = dict(
        confidence_score=c["confidence_score"], ptm=c["ptm"], iptm=c["iptm"]
    )

# top-confidence model, used for the visual overlay
top_model_idx = max(conf_by_model, key=lambda k: conf_by_model[k]["confidence_score"])

results = []


def load_ref():
    cmd.load(REF_PDB, "ref")
    cmd.remove("ref and solvent")
    cmd.remove("ref and hydro")
    cmd.remove("ref and not alt ''+A")
    cmd.alter("ref", "alt=''")
    cmd.remove("ref and not (chain A+D)")
    cmd.alter("ref and chain D", "resi=str(int(resi)+4)")
    cmd.sort("ref")
    cmd.rebuild()


for i in range(N_MODELS):
    cmd.reinitialize()
    load_ref()

    cif = f"{BOLTZ_DIR}/cSrc_SH2_EPQpYEEIPIYL_model_{i}.cif"
    cmd.load(cif, "bz")
    cmd.remove("bz and hydro")

    cmd.create("bz_pepfit", "bz")

    mobile_sel = "bz and chain A and polymer.protein and name N+CA+C+O"
    target_sel = "ref and chain A and polymer.protein and name N+CA+C+O"
    aln = cmd.align(mobile_sel, target_sel, cycles=5, object=f"aln_prot_{i}")
    rmsd_after, natoms_after, ncycles, rmsd_before, natoms_before, score, nres = aln

    resi_str = "+".join(str(r) for r in PEP_RESI)
    bz_pep_bb = f"bz and chain B and resi {resi_str} and name N+CA+C+O"
    ref_pep_bb = f"ref and chain D and resi {resi_str} and name N+CA+C+O"
    n_bz = cmd.count_atoms(bz_pep_bb)
    n_ref = cmd.count_atoms(ref_pep_bb)
    assert n_bz == n_ref, f"model_{i}: backbone atom count mismatch bz={n_bz} ref={n_ref}"
    pep_bb_rmsd_inplace = cmd.rms_cur(bz_pep_bb, ref_pep_bb, matchmaker=-1)

    bz_pep_ca = f"bz and chain B and resi {resi_str} and name CA"
    ref_pep_ca = f"ref and chain D and resi {resi_str} and name CA"
    pep_ca_rmsd_inplace = cmd.rms_cur(bz_pep_ca, ref_pep_ca, matchmaker=-1)

    bz_ptr = f"bz and chain B and resi {PTR_RESI} and name {PTR_SIDECHAIN_ATOMS}"
    ref_ptr = f"ref and chain D and resi {PTR_RESI} and name {PTR_SIDECHAIN_ATOMS}"
    n_bz_ptr = cmd.count_atoms(bz_ptr)
    n_ref_ptr = cmd.count_atoms(ref_ptr)
    if n_bz_ptr == n_ref_ptr and n_bz_ptr > 0:
        ptr_rmsd = cmd.rms_cur(bz_ptr, ref_ptr, matchmaker=-1)
    else:
        ptr_rmsd = float("nan")

    per_res = []
    for resi in range(1, 105):
        bz_ca_sel = f"bz and chain A and resi {resi} and name CA"
        ref_ca_sel = f"ref and chain A and resi {resi} and name CA"
        if cmd.count_atoms(bz_ca_sel) == 1 and cmd.count_atoms(ref_ca_sel) == 1:
            c1 = cmd.get_atom_coords(bz_ca_sel)
            c2 = cmd.get_atom_coords(ref_ca_sel)
            d = float(np.linalg.norm(np.array(c1) - np.array(c2)))
            per_res.append((resi, d))
    dists = np.array([d for _, d in per_res])
    mean_d = float(dists.mean())
    max_d = float(dists.max())
    max_resi = per_res[int(dists.argmax())][0]
    outliers = [(r, d) for r, d in per_res if d > mean_d + 2 * dists.std()]
    outliers_str = ", ".join(
        f"resi {r} ({d:.2f} A)" for r, d in sorted(outliers, key=lambda x: -x[1])[:8]
    )

    bzpf_pep_bb = f"bz_pepfit and chain B and resi {resi_str} and name N+CA+C+O"
    aln_pep = cmd.align(bzpf_pep_bb, ref_pep_bb, cycles=5, object=f"aln_pep_{i}")
    pep_local_rmsd_after, pep_local_n_after, _, pep_local_rmsd_before, pep_local_n_before, _, _ = aln_pep

    row = dict(
        model=f"model_{i}",
        confidence_score=conf_by_model[i]["confidence_score"],
        ptm=conf_by_model[i]["ptm"],
        iptm=conf_by_model[i]["iptm"],
        protein_rmsd_before=rmsd_before,
        protein_natoms_before=natoms_before,
        protein_rmsd_after=rmsd_after,
        protein_natoms_after=natoms_after,
        peptide_bb_rmsd_inplace=pep_bb_rmsd_inplace,
        peptide_ca_rmsd_inplace=pep_ca_rmsd_inplace,
        peptide_local_rmsd_before=pep_local_rmsd_before,
        peptide_local_rmsd_after=pep_local_rmsd_after,
        ptr_sidechain_rmsd=ptr_rmsd,
        per_res_mean=mean_d,
        per_res_max=max_d,
        per_res_max_resi=max_resi,
        outliers=outliers_str,
    )
    results.append(row)
    print(
        f"[model_{i:2d}] conf={row['confidence_score']:.4f} "
        f"prot_bb={rmsd_after:.3f}A pep_bb_inplace={pep_bb_rmsd_inplace:.3f}A "
        f"pep_local={pep_local_rmsd_after:.3f}A pTyr_sc={ptr_rmsd:.3f}A "
        f"CA_max={max_d:.2f}A@resi{max_resi}"
    )

    if i == top_model_idx:
        cmd.hide("everything")
        cmd.bg_color("white")
        cmd.show("cartoon", "bz and chain A or ref and chain A")
        cmd.color("palegreen", "bz and chain A")
        cmd.color("salmon", "ref and chain A")
        cmd.show("sticks", "bz and chain B and not hydro")
        cmd.show("sticks", "ref and chain D and not hydro")
        cmd.color("green", "bz and chain B and elem C")
        cmd.color("red", "ref and chain D and elem C")
        cmd.set("cartoon_transparency", 0.2)
        cmd.util.cnc("bz and chain B")
        cmd.util.cnc("ref and chain D")
        cmd.set("ray_opaque_background", 0)
        cmd.orient("bz and chain A or ref and chain A")
        cmd.zoom("bz and chain A or ref and chain A or bz and chain B or ref and chain D", buffer=3)
        cmd.set("ray_trace_mode", 1)
        cmd.ray(1600, 1600)
        cmd.png(OUT_PNG, dpi=300)
        cmd.save(OUT_PSE)
        print(f"[model_{i}] (top confidence) saved overlay -> {OUT_PNG} / {OUT_PSE}")

    cmd.delete("bz_pepfit")

print("\n\n===================== SUMMARY (all 25 models) =====================")
header = (
    f"{'model':10s} {'conf':>7s} {'prot_bb(A)':>11s} "
    f"{'pep_bb_inplace(A)':>18s} {'pep_CA_inplace(A)':>18s} "
    f"{'pep_local(A)':>13s} {'pTyr_sc(A)':>11s} {'CA_max(A)':>10s}"
)
print(header)
for r in results:
    print(
        f"{r['model']:10s} {r['confidence_score']:7.4f} {r['protein_rmsd_after']:11.3f} "
        f"{r['peptide_bb_rmsd_inplace']:18.3f} {r['peptide_ca_rmsd_inplace']:18.3f} "
        f"{r['peptide_local_rmsd_after']:13.3f} {r['ptr_sidechain_rmsd']:11.3f} "
        f"{r['per_res_max']:10.3f}"
    )

arr = lambda key: np.array([r[key] for r in results])
print("\n---- ensemble stats across all 25 samples (mean +/- std [min, max]) ----")
for key, label in [
    ("protein_rmsd_after", "protein backbone RMSD"),
    ("peptide_bb_rmsd_inplace", "peptide backbone RMSD in-place"),
    ("peptide_ca_rmsd_inplace", "peptide CA RMSD in-place"),
    ("peptide_local_rmsd_after", "peptide local (self-fit) RMSD"),
    ("ptr_sidechain_rmsd", "pTyr side-chain RMSD"),
]:
    a = arr(key)
    print(f"{label:35s}: {a.mean():.3f} +/- {a.std():.3f}  [{a.min():.3f}, {a.max():.3f}]")

with open(OUT_JSON, "w") as fh:
    json.dump(dict(top_model_idx=top_model_idx, results=results), fh, indent=2)
print(f"\nFull results JSON -> {OUT_JSON}")

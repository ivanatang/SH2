"""
align_af3_vs_1sps_5models.py

Structural validation: how well do the 5 AlphaFold3 models of WT c-Src SH2 +
phosphopeptide pYEEI (EPQ-pY-EEIPIYL) reproduce the experimental crystal
structure 1SPS (Waksman & Kuriyan 1993, 2.7 A, Src SH2 + pYEEI peptide)?

Inputs
------
- AF3 models (mmCIF), chain A = SH2 protein (104 res), chain B = peptide
  (11 res, PTR at resi 4):
    /Users/ivanatang/Developer/SH2/af3/fold_wt_sh2_pyeei/fold_wt_sh2_pyeei_model_{0..4}.cif
- Crystal reference:
    /Users/ivanatang/Developer/SH2/structures/1SPS.pdb
  Chains A/B/C = SH2 domain (3 crystallographically independent copies),
  chains D/E/F = the bound pYEEI peptide (3 independent copies). Chain A/D is
  used as the primary reference pair per the task spec.

IMPORTANT numbering note (verified by direct inspection of both files)
------------------------------------------------------------------------
1SPS chain D (the primary reference peptide copy) only has *7 of the 11*
peptide residues resolved in electron density: PDB-numbered PRO(-2), GLN(-1),
PTR(0), GLU(1), GLU(2), ILE(3), PRO(4) -- i.e. SEQRES/peptide positions 2-8.
The N-terminal Glu (position 1) and the C-terminal Ile-Tyr-Leu (positions
9-11) are disordered/not modeled in chain D. (Chain E is resolved positions
2-7 only -- even less complete; chain F is the most complete copy, resolved
positions 2-11, missing only the N-terminal Glu.) This script therefore
restricts all peptide-vs-peptide RMSD comparisons to the AF3 peptide residues
2-8, which is the only region that can be directly compared against the
chain-D reference. This is itself a notable result worth reporting: AF3
predicts full-length ordered density for a region that is crystallographically
disordered in 2/3 copies of the ASU, so agreement in the ordered "core"
(pTyr pocket + immediate flanks) is the meaningful test; positions 1 and 9-11
cannot be judged against chain D at all (see chain F for a partial check on
position 9-11 if desired -- not automated here).

To directly compare atom-for-atom, 1SPS chain D peptide residues are
renumbered in-place (+4) onto the AF3 peptide's 1-11 numbering frame before
any RMSD calculations (chain D PDB resi -2..4 -> 2..8, matching AF3 chain B
resi 2..8, with PTR landing on resi 4 in both).

Protein chain A numbering is already identical (1-104) between AF3 and 1SPS
(SEQRES-verified identical sequence), so no protein renumbering is needed.

Atom order within each residue (N, CA, C, O, ...) was verified identical and
ascending-residue in both files, so positional atom pairing (matchmaker=-1)
is valid for cmd.rms_cur() on the backbone and PTR side-chain selections
below; atom *names* for the PTR side chain (CB,CG,CD1,CD2,CE1,CE2,CZ,OH,P,
O1P,O2P,O3P) were also verified identical between AF3 and 1SPS.

What this script does, per model (5x)
--------------------------------------
1. Loads 1SPS, keeps only chain A (protein) + chain D (peptide), strips
   solvent/altlocs, renumbers chain D peptide onto AF3's 1-11 frame.
2. Loads the AF3 model, superposes ONLY the SH2 protein backbone (N,CA,C,O)
   onto 1SPS chain A via cmd.align (5 cycles, outlier rejection) -- this
   transforms the whole AF3 object (protein + peptide together).
3. Reports protein backbone RMSD before/after outlier rejection.
4. Without further fitting, computes peptide backbone RMSD "in place"
   (cmd.rms_cur, no refit) for both N,CA,C,O and CA-only, over the resolved
   region (AF3 resi 2-8) -- this is the key docking-pose accuracy metric.
5. Independently, superposes ONLY the peptide backbone (AF3 peptide onto
   1SPS chain D peptide, on a duplicated/untransformed copy of the AF3
   object) to isolate local peptide conformational accuracy from global
   docking pose.
6. Reports PTR side-chain heavy-atom RMSD after the protein-based
   superposition (step 2's transform, no refit).
7. Per-residue CA deviation on the SH2 protein (post protein-fit) to flag
   any loop regions (e.g. BC loop, EF loop) with large local disagreement.
8. Saves an overlay PNG/PSE for the top-ranked model_0.

Usage
-----
    pymol -cq /Users/ivanatang/Developer/SH2/af3/align_af3_vs_1sps_5models.py

Outputs
-------
    /Users/ivanatang/Developer/SH2/structures/af3_vs_1sps_model0_overlay.png
    /Users/ivanatang/Developer/SH2/structures/af3_vs_1sps_model0_overlay.pse
    /Users/ivanatang/Developer/SH2/structures/af3_vs_1sps_results.json
plus a full numeric summary table printed to stdout/log.
"""

import json

import numpy as np
from pymol import cmd

BASE = "/Users/ivanatang/Developer/SH2"
REF_PDB = f"{BASE}/structures/1SPS.pdb"
AF3_DIR = f"{BASE}/af3/fold_wt_sh2_pyeei"
OUT_PNG = f"{BASE}/structures/af3_vs_1sps_model0_overlay.png"
OUT_PSE = f"{BASE}/structures/af3_vs_1sps_model0_overlay.pse"
OUT_JSON = f"{BASE}/structures/af3_vs_1sps_results.json"

N_MODELS = 5

# AF3 peptide: chain B, resi 1-11 (EPQ-[PTR]-EEIPIYL), PTR at resi 4.
# 1SPS chain D (primary ref copy) is only resolved for peptide positions
# 2-8 (see numbering note above); after renumbering chain D by +4 those
# same positions read as resi 2-8 in both objects.
PEP_RESI = list(range(2, 9))  # 2..8, the resolved/comparable region
PTR_RESI = 4

PTR_SIDECHAIN_ATOMS = "CB+CG+CD1+CD2+CE1+CE2+CZ+OH+P+O1P+O2P+O3P"

results = []


def load_ref():
    cmd.load(REF_PDB, "ref")
    cmd.remove("ref and solvent")
    cmd.remove("ref and hydro")
    cmd.remove("ref and not alt ''+A")
    cmd.alter("ref", "alt=''")
    cmd.remove("ref and not (chain A+D)")
    # renumber peptide chain D onto AF3's 1-11 numbering frame
    cmd.alter("ref and chain D", "resi=str(int(resi)+4)")
    cmd.sort("ref")
    cmd.rebuild()


for i in range(N_MODELS):
    cmd.reinitialize()
    load_ref()

    cif = f"{AF3_DIR}/fold_wt_sh2_pyeei_model_{i}.cif"
    cmd.load(cif, "af3")
    cmd.remove("af3 and hydro")

    af3_chains = set()
    cmd.iterate("af3", "af3_chains.add(chain)", space={"af3_chains": af3_chains})
    ref_chains = set()
    cmd.iterate("ref", "ref_chains.add(chain)", space={"ref_chains": ref_chains})
    print(f"[model_{i}] af3 chains={sorted(af3_chains)}  ref chains={sorted(ref_chains)}")

    # duplicate af3 in its ORIGINAL (untransformed) frame, used later for the
    # independent local peptide-only superposition (step 5)
    cmd.create("af3_pepfit", "af3")

    # ---------------- STEP 2/3: global protein-backbone superposition ----------------
    mobile_sel = "af3 and chain A and polymer.protein and name N+CA+C+O"
    target_sel = "ref and chain A and polymer.protein and name N+CA+C+O"
    aln = cmd.align(mobile_sel, target_sel, cycles=5, object=f"aln_prot_{i}")
    rmsd_after, natoms_after, ncycles, rmsd_before, natoms_before, score, nres = aln
    print(
        f"[model_{i}] protein align: before={rmsd_before:.3f} A ({natoms_before} at), "
        f"after={rmsd_after:.3f} A ({natoms_after} at), cycles={ncycles}, nres_aligned={nres}"
    )

    # ---------------- STEP 4: peptide backbone RMSD in place (no refit) ----------------
    resi_str = "+".join(str(r) for r in PEP_RESI)
    af3_pep_bb = f"af3 and chain B and resi {resi_str} and name N+CA+C+O"
    ref_pep_bb = f"ref and chain D and resi {resi_str} and name N+CA+C+O"
    n_af3 = cmd.count_atoms(af3_pep_bb)
    n_ref = cmd.count_atoms(ref_pep_bb)
    assert n_af3 == n_ref, f"backbone atom count mismatch: af3={n_af3} ref={n_ref}"
    pep_bb_rmsd_inplace = cmd.rms_cur(af3_pep_bb, ref_pep_bb, matchmaker=-1)

    af3_pep_ca = f"af3 and chain B and resi {resi_str} and name CA"
    ref_pep_ca = f"ref and chain D and resi {resi_str} and name CA"
    pep_ca_rmsd_inplace = cmd.rms_cur(af3_pep_ca, ref_pep_ca, matchmaker=-1)

    print(
        f"[model_{i}] peptide backbone RMSD in-place (N,CA,C,O; resi {PEP_RESI[0]}-{PEP_RESI[-1]}): "
        f"{pep_bb_rmsd_inplace:.3f} A ({n_af3} atoms);  CA-only: {pep_ca_rmsd_inplace:.3f} A"
    )

    # ---------------- STEP 6: pTyr side-chain heavy-atom RMSD (no refit) ----------------
    af3_ptr = f"af3 and chain B and resi {PTR_RESI} and name {PTR_SIDECHAIN_ATOMS}"
    ref_ptr = f"ref and chain D and resi {PTR_RESI} and name {PTR_SIDECHAIN_ATOMS}"
    n_af3_ptr = cmd.count_atoms(af3_ptr)
    n_ref_ptr = cmd.count_atoms(ref_ptr)
    if n_af3_ptr == n_ref_ptr and n_af3_ptr > 0:
        ptr_rmsd = cmd.rms_cur(af3_ptr, ref_ptr, matchmaker=-1)
    else:
        print(f"[model_{i}] WARNING: PTR sidechain atom count mismatch af3={n_af3_ptr} ref={n_ref_ptr}")
        ptr_rmsd = float("nan")
    print(f"[model_{i}] pTyr side-chain heavy-atom RMSD (no refit): {ptr_rmsd:.3f} A ({n_af3_ptr} atoms)")

    # ---------------- STEP 7: per-residue CA deviation on the SH2 protein ----------------
    per_res = []
    for resi in range(1, 105):
        af3_ca_sel = f"af3 and chain A and resi {resi} and name CA"
        ref_ca_sel = f"ref and chain A and resi {resi} and name CA"
        if cmd.count_atoms(af3_ca_sel) == 1 and cmd.count_atoms(ref_ca_sel) == 1:
            c1 = cmd.get_atom_coords(af3_ca_sel)
            c2 = cmd.get_atom_coords(ref_ca_sel)
            d = float(np.linalg.norm(np.array(c1) - np.array(c2)))
            per_res.append((resi, d))
    dists = np.array([d for _, d in per_res])
    mean_d = float(dists.mean())
    max_d = float(dists.max())
    max_resi = per_res[int(dists.argmax())][0]
    outliers = [(r, d) for r, d in per_res if d > mean_d + 2 * dists.std()]
    outliers_str = ", ".join(f"resi {r} ({d:.2f} A)" for r, d in sorted(outliers, key=lambda x: -x[1])[:8])
    print(f"[model_{i}] per-residue CA (post protein-fit): mean={mean_d:.2f} A, max={max_d:.2f} A at resi {max_resi}")
    print(f"[model_{i}] outliers (>mean+2sd): {outliers_str if outliers_str else 'none'}")

    # ---------------- STEP 5: independent LOCAL peptide-only superposition ----------------
    af3pf_pep_bb = f"af3_pepfit and chain B and resi {resi_str} and name N+CA+C+O"
    aln_pep = cmd.align(af3pf_pep_bb, ref_pep_bb, cycles=5, object=f"aln_pep_{i}")
    pep_local_rmsd_after, pep_local_n_after, _, pep_local_rmsd_before, pep_local_n_before, _, _ = aln_pep
    print(
        f"[model_{i}] peptide-only LOCAL superposition: before={pep_local_rmsd_before:.3f} A "
        f"({pep_local_n_before} at), after={pep_local_rmsd_after:.3f} A ({pep_local_n_after} at)"
    )

    results.append(
        dict(
            model=f"model_{i}",
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
    )

    # ---------------- visual sanity check for top-ranked model_0 ----------------
    if i == 0:
        cmd.hide("everything")
        cmd.bg_color("white")
        cmd.show("cartoon", "af3 and chain A or ref and chain A")
        cmd.color("skyblue", "af3 and chain A")
        cmd.color("salmon", "ref and chain A")
        cmd.show("sticks", "af3 and chain B and not hydro")
        cmd.show("sticks", "ref and chain D and not hydro")
        cmd.color("blue", "af3 and chain B and elem C")
        cmd.color("red", "ref and chain D and elem C")
        cmd.set("cartoon_transparency", 0.2)
        cmd.util.cnc("af3 and chain B")
        cmd.util.cnc("ref and chain D")
        cmd.set("ray_opaque_background", 0)
        cmd.orient("af3 and chain A or ref and chain A")
        cmd.zoom("af3 and chain A or ref and chain A or af3 and chain B or ref and chain D", buffer=3)
        cmd.set("ray_trace_mode", 1)
        cmd.ray(1600, 1600)
        cmd.png(OUT_PNG, dpi=300)
        cmd.save(OUT_PSE)
        print(f"[model_0] saved overlay image -> {OUT_PNG}")
        print(f"[model_0] saved session      -> {OUT_PSE}")

    cmd.delete("af3_pepfit")

# ---------------------------- summary table ----------------------------
print("\n\n===================== SUMMARY =====================")
header = (
    f"{'model':10s} {'prot_bb_RMSD(A)':>16s} {'prot_Natoms':>12s} "
    f"{'pep_bb_inplace(A)':>18s} {'pep_CA_inplace(A)':>18s} "
    f"{'pep_local(A)':>13s} {'pTyr_sc(A)':>11s} {'CA_max(A)':>10s}"
)
print(header)
for r in results:
    print(
        f"{r['model']:10s} {r['protein_rmsd_after']:16.3f} {r['protein_natoms_after']:12d} "
        f"{r['peptide_bb_rmsd_inplace']:18.3f} {r['peptide_ca_rmsd_inplace']:18.3f} "
        f"{r['peptide_local_rmsd_after']:13.3f} {r['ptr_sidechain_rmsd']:11.3f} "
        f"{r['per_res_max']:10.3f}"
    )

with open(OUT_JSON, "w") as fh:
    json.dump(results, fh, indent=2)
print(f"\nFull results JSON -> {OUT_JSON}")

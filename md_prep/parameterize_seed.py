"""
Parameterize + solvate one AF3/Boltz-2 superbinder-SH2 + pYEEI seed structure
for GROMACS, following exactly the same methodology as
SH2_superbinder_pTyr_peptide_parameterization.ipynb (dianionic phosphotyrosine,
capped ACE/NME peptide termini, ff14SB protein + Sage peptide force-field
split via OpenFF Interchange, NAGL partial charges, rhombic dodecahedron
solvent box) -- generalized to run once per seed instead of interactively for
one system.

Generalization vs. the notebook (important): the notebook was written for a
different peptide (DHEPIYEQWGW, pTyr at residue 6) than the pYEEI peptide our
8 seeds use (EPQYEEIPIYL, pTyr at residue 4), and hardcoded residue numbers
1/6/11 for the termini/phosphate. This script auto-detects the PTR residue
and N-/C-termini from the input structure instead, so it isn't silently
wrong about which atom is the phosphate.

Usage:
    /Users/ivanatang/miniforge3/envs/SH2/bin/python md_prep/parameterize_seed.py \
        --input-cif <path/to/seed.cif> --seed-name af3_c1 --outdir md_prep/seeds/af3_c1

Run once per seed (8x for this project: af3_c1-4, boltz_c1-4).
"""
import argparse
import os
import subprocess
import sys

import numpy as np
from openff.interchange import Interchange
from openff.interchange.components._packmol import RHOMBIC_DODECAHEDRON, pack_box
from openff.toolkit import ForceField, Molecule, Topology
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from openff.units import unit
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

# same peptide sequence (EPQYEEIPIYL: Glu1-Pro2-Gln3-pTyr4-Glu5-Glu6-Ile7-Pro8-
# Ile9-Tyr10-Leu11) across all 8 seeds -- 3 Glu (-1 each) + pTyr dianion (-2),
# no Asp/Lys/Arg/His-relevant charges at pH 7.5. Used only as a sanity check,
# not to hardcode any atom selection.
EXPECTED_PEPTIDE_CHARGE = -5


def find_atom(mol, resnum, name):
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info and info.GetResidueNumber() == resnum and info.GetName().strip() == name:
            return atom.GetIdx()
    raise ValueError(f"atom resnum={resnum} name={name} not found")


def find_ptr_resnum(mol):
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info and info.GetResidueName().strip() == "PTR":
            return info.GetResidueNumber()
    raise ValueError("no PTR residue found in peptide")


def get_resnum_range(mol):
    resnums = set()
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info:
            resnums.add(info.GetResidueNumber())
    return min(resnums), max(resnums)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-cif", required=True)
    ap.add_argument("--seed-name", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    input_cif = os.path.abspath(args.input_cif)

    chainA_pdb = f"{outdir}/{args.seed_name}_chainA_protein_raw.pdb"
    chainB_pdb = f"{outdir}/{args.seed_name}_chainB_peptide_raw.pdb"

    print(f"[{args.seed_name}] splitting {input_cif} into protein/peptide PDBs")
    subprocess.run(
        [
            "pymol", "-cq", "-d",
            f"load {input_cif}, af3; "
            f"save {chainA_pdb}, af3 and chain A and polymer.protein; "
            f"save {chainB_pdb}, af3 and chain B",
        ],
        check=True,
    )

    # ---------------- peptide: dianionic PTR, ACE/NME caps ----------------
    print(f"[{args.seed_name}] building peptide (PTR dianion + ACE/NME caps)")
    peptide_raw = Chem.MolFromPDBBlock(open(chainB_pdb).read(), removeHs=False, sanitize=False)
    Chem.SanitizeMol(peptide_raw)

    n_term_resnum, c_term_resnum = get_resnum_range(peptide_raw)
    ptr_resnum = find_ptr_resnum(peptide_raw)
    print(f"[{args.seed_name}] peptide spans residues {n_term_resnum}-{c_term_resnum}, PTR at {ptr_resnum}")

    rw = Chem.RWMol(peptide_raw)
    p_idx, o1p_idx, o2p_idx, o3p_idx = (find_atom(rw, ptr_resnum, n) for n in ("P", "O1P", "O2P", "O3P"))
    for idx in (p_idx, o1p_idx, o2p_idx, o3p_idx):
        rw.GetAtomWithIdx(idx).SetFormalCharge(0)
        rw.GetAtomWithIdx(idx).SetNoImplicit(False)
    rw.GetBondBetweenAtoms(p_idx, o3p_idx).SetBondType(Chem.BondType.DOUBLE)
    rw.GetAtomWithIdx(o1p_idx).SetFormalCharge(-1)
    rw.GetAtomWithIdx(o2p_idx).SetFormalCharge(-1)
    for idx in (p_idx, o1p_idx, o2p_idx, o3p_idx):
        rw.GetAtomWithIdx(idx).SetNoImplicit(True)

    try:
        rw.RemoveAtom(find_atom(rw, c_term_resnum, "OXT"))
    except ValueError:
        pass

    n_fixed = rw.GetNumAtoms()

    n_idx = find_atom(rw, n_term_resnum, "N")
    ace_c, ace_o, ace_ch3 = (rw.AddAtom(Chem.Atom(e)) for e in ("C", "O", "C"))
    rw.AddBond(n_idx, ace_c, Chem.BondType.SINGLE)
    rw.AddBond(ace_c, ace_o, Chem.BondType.DOUBLE)
    rw.AddBond(ace_c, ace_ch3, Chem.BondType.SINGLE)
    rw.GetAtomWithIdx(n_idx).SetNoImplicit(True)
    rw.GetAtomWithIdx(n_idx).SetNumExplicitHs(1)

    c_idx = find_atom(rw, c_term_resnum, "C")
    nme_n, nme_ch3 = (rw.AddAtom(Chem.Atom(e)) for e in ("N", "C"))
    rw.AddBond(c_idx, nme_n, Chem.BondType.SINGLE)
    rw.AddBond(nme_n, nme_ch3, Chem.BondType.SINGLE)
    rw.GetAtomWithIdx(nme_n).SetNoImplicit(True)
    rw.GetAtomWithIdx(nme_n).SetNumExplicitHs(1)

    peptide = rw.GetMol()
    Chem.SanitizeMol(peptide)
    charge = Chem.GetFormalCharge(peptide)
    print(f"[{args.seed_name}] peptide net formal charge: {charge} (expected {EXPECTED_PEPTIDE_CHARGE})")
    assert charge == EXPECTED_PEPTIDE_CHARGE, "unexpected net charge -- check residue protonation assumptions"

    coord_map = {i: peptide.GetConformer().GetAtomPosition(i) for i in range(n_fixed)}
    # the notebook this is based on notes that cap embedding can fail and
    # suggests trying a different randomSeed -- that alone isn't reliable
    # (verified: 50 seeds x 200 attempts each still failed for one cluster's
    # peptide). Root cause: default ETKDG adds idealized bond-length/angle/
    # planarity "chemical knowledge" bounds on top of pure distance geometry,
    # which can be inconsistent with a coordMap of REAL (non-idealized) AF3
    # coordinates -- making the combined bounds matrix infeasible even though
    # the underlying molecular graph embeds fine on its own (verified directly).
    # Disabling useBasicKnowledge relaxes those idealized bounds so the
    # pinned real coordinates stay feasible; try random seeds within that
    # relaxed setting first, only falling back to the notebook's original
    # (stricter) settings if that somehow doesn't work either.
    from rdkit.Chem import rdDistGeom

    cid = -1
    for seed in (42, 7, 123, 2024, 31415):
        params = rdDistGeom.ETKDGv3()
        params.useRandomCoords = True
        params.randomSeed = seed
        params.SetCoordMap(coord_map)
        params.enforceChirality = False
        params.useBasicKnowledge = False
        params.maxIterations = 2000
        cid = rdDistGeom.EmbedMolecule(peptide, params)
        if cid != -1:
            print(f"[{args.seed_name}] cap embedding succeeded (relaxed ETKDGv3) with randomSeed={seed}")
            break
    if cid == -1:
        for seed in (42, 7, 123, 2024, 31415):
            cid = AllChem.EmbedMolecule(peptide, coordMap=coord_map, useRandomCoords=True, randomSeed=seed, maxAttempts=200)
            if cid != -1:
                print(f"[{args.seed_name}] cap embedding succeeded (default ETKDG) with randomSeed={seed}")
                break
    assert cid != -1, "cap embedding failed with all tried random seeds and both ETKDG configurations"

    coords = peptide.GetConformer().GetPositions()
    max_drift = max(
        np.linalg.norm(coords[i] - np.array([coord_map[i].x, coord_map[i].y, coord_map[i].z]))
        for i in range(n_fixed)
    )
    assert max_drift < 1e-3, "embedding moved atoms it should have held fixed"
    min_dist = min(
        np.linalg.norm(coords[i] - coords[j])
        for i in range(len(coords)) for j in range(i + 1, len(coords))
        if not peptide.GetBondBetweenAtoms(i, j)
    )
    assert min_dist > 1.5, "cap placement clashes with an existing atom"

    peptide = Chem.AddHs(peptide, addCoords=True)
    peptide_offmol = Molecule.from_rdkit(peptide, allow_undefined_stereo=True)
    peptide_offmol.name = "pTyr_peptide"
    for atom in peptide_offmol.atoms:
        atom.metadata["residue_name"] = "LIG"
        atom.metadata["residue_number"] = "1"

    print(f"[{args.seed_name}] assigning NAGL partial charges")
    peptide_offmol.assign_partial_charges(
        partial_charge_method="openff-gnn-am1bcc-1.0.0.pt",
        toolkit_registry=NAGLToolkitWrapper(),
    )
    charge_sum = sum(q.m for q in peptide_offmol.partial_charges)
    assert abs(charge_sum - peptide_offmol.total_charge.m) < 1e-3

    peptide_intrcg = Interchange.from_smirnoff(
        force_field=ForceField("openff_unconstrained-2.0.0.offxml"),
        topology=[peptide_offmol],
        charge_from_molecules=[peptide_offmol],
    )

    # ---------------- protein: pdbfixer + ff14SB ----------------
    print(f"[{args.seed_name}] fixing protein (missing atoms/H at pH 7.5) + ff14SB")
    fixer = PDBFixer(filename=chainA_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.5)
    protein_fixed_pdb = f"{outdir}/{args.seed_name}_chainA_protein_fixed_H.pdb"
    with open(protein_fixed_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    protein_full = Topology.from_pdb(protein_fixed_pdb)
    protein = protein_full.molecule(0)
    protein.name = "protein"

    ff14sb = ForceField("ff14sb_off_impropers_0.0.3.offxml")
    protein_intrcg = Interchange.from_smirnoff(force_field=ff14sb, topology=protein.to_topology())

    # ---------------- dock + solvate ----------------
    print(f"[{args.seed_name}] docking protein+peptide, building solvent box")
    docked_intrcg = protein_intrcg.combine(peptide_intrcg)
    total_charge = round(sum(docked_intrcg["Electrostatics"].charges.values()), 3)
    assert total_charge == protein.total_charge + peptide_offmol.total_charge

    total_charge_e = float(total_charge.m)
    if total_charge_e < 0:
        counterion = Molecule.from_smiles("[Na+]")
        counterion.name = "NA"
        n_counterions = int(round(abs(total_charge_e)))
    elif total_charge_e > 0:
        counterion = Molecule.from_smiles("[Cl-]")
        counterion.name = "CL"
        n_counterions = int(round(abs(total_charge_e)))
    else:
        counterion = None
        n_counterions = 0

    water = Molecule.from_smiles("O")
    water.name = "SOL"
    water.generate_conformers(n_conformers=1)
    print(f"[{args.seed_name}] net charge: {total_charge_e}, counterion: {counterion.name if counterion else None} x {n_counterions}")

    xyz = protein.conformers[0].to(unit.nanometer).m
    centroid = xyz.mean(axis=0)
    protein_radius_nm = np.sqrt(((xyz - centroid) ** 2).sum(axis=1).max())
    buffer_nm = 2.0
    scale_nm = 2.0 * protein_radius_nm + buffer_nm
    box_vectors = (scale_nm * RHOMBIC_DODECAHEDRON) * unit.nanometer

    box_nm = box_vectors.to(unit.nanometer).m
    V_box_nm3 = abs(np.linalg.det(box_nm))
    V_solute_nm3 = (4.0 / 3.0) * np.pi * (protein_radius_nm ** 3)
    waters_per_nm3 = 33.4
    n_water = int(waters_per_nm3 * max(V_box_nm3 - V_solute_nm3, 0.0))
    print(f"[{args.seed_name}] box scale {scale_nm:.3f} nm, packing {n_water} waters, {n_counterions} ions")

    molecules = [water]
    number_of_copies = [n_water]
    if counterion is not None and n_counterions > 0:
        molecules.append(counterion)
        number_of_copies.append(n_counterions)

    packmol_dir = f"{outdir}/packmol_solv_dodecahedron"
    packed_topology = pack_box(
        solute=docked_intrcg.topology,
        molecules=molecules,
        number_of_copies=number_of_copies,
        box_vectors=box_vectors,
        center_solute=True,
        tolerance=2.0 * unit.angstrom,
        working_directory=packmol_dir,
        retain_working_files=False,
    )
    print(f"[{args.seed_name}] packed {packed_topology.n_molecules} molecules")

    topology_molecules = [water] * n_water
    if counterion is not None:
        topology_molecules += [counterion] * n_counterions
    water_intrcg = Interchange.from_smirnoff(
        force_field=ForceField("openff_unconstrained-2.0.0.offxml"),
        topology=topology_molecules,
    )

    system_intrcg = docked_intrcg.combine(water_intrcg)
    system_intrcg.positions = packed_topology.get_positions()
    system_intrcg.box = packed_topology.box_vectors

    gromacs_dir = f"{outdir}/gromacs"
    os.makedirs(gromacs_dir, exist_ok=True)
    prefix = f"{args.seed_name}_dodecahedron"
    cwd = os.getcwd()
    os.chdir(gromacs_dir)
    try:
        system_intrcg.to_gromacs(prefix=prefix, decimal=3, monolithic=False)

        # Interchange's GROMACS export doesn't generate backbone position
        # restraints (used by NVT/NPT via `define = -DPOSRES`) -- the reference
        # system's posre_protein.itp was added as a separate manual step, so
        # reproduce that here: restrain the Backbone group (N/CA/C), not all
        # protein atoms, matching the reference file exactly.
        print(f"[{args.seed_name}] generating posre_protein.itp (Backbone group)")
        gro = f"{prefix}.gro"
        subprocess.run(["gmx", "make_ndx", "-f", gro, "-o", "index.ndx"], input="q\n", text=True, capture_output=True, check=True)
        subprocess.run(
            ["gmx", "genrestr", "-f", gro, "-n", "index.ndx", "-o", "posre_protein.itp"],
            input="Backbone\n", text=True, capture_output=True, check=True,
        )
        os.remove("index.ndx")

        top_path = f"{prefix}.top"
        with open(top_path) as fh:
            top = fh.read()
        marker = f'#include "{prefix}_protein.itp"\n'
        assert marker in top, "protein.itp include not found in .top -- can't insert POSRES block"
        top = top.replace(marker, marker + '#ifdef POSRES\n#include "posre_protein.itp"\n#endif\n')
        with open(top_path, "w") as fh:
            fh.write(top)
    finally:
        os.chdir(cwd)

    print(f"[{args.seed_name}] DONE -> {gromacs_dir}/{prefix}.gro / .top")


if __name__ == "__main__":
    main()

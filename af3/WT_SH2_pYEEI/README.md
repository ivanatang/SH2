# AF3 inputs

`superbinder_sh2_pTyr.json` — chain A sequence is the superbinder SH2 domain (UniProt P12931
numbering, residues 144–252, T183V/C188A/K206L), sourced from `structures/4F5B_superbinder_SH2_pTyr.pdb`
chain A (verified 2026-08-03: extracted directly from ATOM records, byte-identical across all three
deposited superbinder structures — 4F59, 4F5A, 4F5B — since they're the same construct with different
bound ligands). Residue 252 (Lys) is present in SEQRES but unresolved in all three crystal structures;
included here from SEQRES since it's part of the actual construct.

Chain B is a placeholder peptide (c-Src consensus `DHEPIYEQWGW` with PTR at the central Tyr) — swap for
the actual peptide of interest.

## Output numbering gotcha

The JSON input above gives chain A as a raw sequence with no residue numbers, so **AF3's output
renumbers chain A starting at 1**, not from 144 (UniProt/crystal numbering used everywhere else in
this project, e.g. in `structures/`). To map between them: `AF3_resi = UniProt_resi - 143`.

Confirmed mutation positions in the AF3 output (verified 2026-08-04, `superbinder_sh2_ptyr_peptide_model.cif`):

| UniProt / crystal resi | AF3 output resi | Residue |
|---|---|---|
| 183 (T183V) | 40 | VAL |
| 188 (C188A) | 45 | ALA |
| 206 (K206L) | 63 | LEU |

All three superbinder mutations are present and correctly predicted. Use the offset above for any
future direct residue lookups against the AF3 `.cif` — selections using crystal-numbering residue
IDs against the AF3 model will silently return 0 atoms rather than erroring.

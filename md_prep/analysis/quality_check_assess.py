#!/usr/bin/env python3
"""
Stage 2 of the production-trajectory quality check (companion: quality_check_extract.sh).

Reads what quality_check_extract.sh wrote under <seed>/prod_md_48x1/qc/ for each
seed and applies the actual pass/warn/fail thresholds for:
  1. atom clashes / blow-ups   (log warnings + potential/temperature trace)
  2. RMSD trace shape          (plateaued fluctuation vs unbounded drift)
  3. peptide dissociation      (protein-peptide min distance + contact count)

All thresholds are domain judgment calls, not physical constants -- they're named
constants below so they're easy to see and adjust, not buried in the logic.

Usage: python3 quality_check_assess.py [--root SCRATCH_ROOT] [seed ...]
  (no seeds = all 8; default root = /scratch/alpine/ivta1597/SH2/seeds)
"""
import argparse
import sys
from pathlib import Path

import numpy as np

ALL_SEEDS = ["af3_c1", "af3_c2", "af3_c3", "af3_c4", "boltz_c1", "boltz_c2", "boltz_c3", "boltz_c4"]
DEFAULT_ROOT = "/scratch/alpine/ivta1597/SH2/seeds"

# --- thresholds (adjust here, not in the logic below) ---
LINCS_WARN_COUNT_FOR_WARN = 5        # a handful of transient LINCS warnings -> WARN
LINCS_WARN_COUNT_FOR_FAIL = 50       # sustained LINCS warnings -> genuine instability -> FAIL
TEMP_MAX_K_FOR_FAIL = 600.0          # target is ~300K; this high = integration blow-up
ENERGY_SPIKE_MAD_MULTIPLE = 15.0     # frame-to-frame potential-energy jump vs typical noise

RMSD_LAST_FRAC = 0.10                # "final" window = last 10% of the trajectory
RMSD_WARN_ANGSTROM = 5.0             # final-window mean backbone RMSD
RMSD_FAIL_ANGSTROM = 8.0
RMSD_SECOND_HALF_SLOPE_WARN = 0.01   # Angstrom/ns, fit over the second half of the run
RMSD_SECOND_HALF_SLOPE_FAIL = 0.02   # still climbing at this rate near the end = not plateaued

DISSOC_LAST_FRAC = 0.10
DISSOC_MINDIST_WARN_NM = 0.6         # min heavy-atom distance, protein vs peptide
DISSOC_MINDIST_FAIL_NM = 0.8
DISSOC_CONTACT_CUTOFF_NM = 0.5       # must match -d in quality_check_extract.sh


def read_xvg(path):
    """Read a GROMACS -xvg none .xvg file (plain columns of numbers) into an (N, ncol) array."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            rows.append([float(x) for x in line.split()])
    return np.array(rows)


def worst(*verdicts):
    order = {"PASS": 0, "WARN": 1, "FAIL": 2}
    return max(verdicts, key=lambda v: order[v])


def check_clashes(qc_dir):
    notes = []
    verdict = "PASS"

    warn_file = qc_dir / "log_warnings.txt"
    lines = warn_file.read_text().splitlines() if warn_file.exists() else []
    n_warnings = len(lines)
    has_fatal = any("Fatal error" in l for l in lines)

    if has_fatal:
        verdict = worst(verdict, "FAIL")
        notes.append("log contains 'Fatal error'")
    elif n_warnings >= LINCS_WARN_COUNT_FOR_FAIL:
        verdict = worst(verdict, "FAIL")
        notes.append(f"{n_warnings} instability warnings in log (>= {LINCS_WARN_COUNT_FOR_FAIL})")
    elif n_warnings >= LINCS_WARN_COUNT_FOR_WARN:
        verdict = worst(verdict, "WARN")
        notes.append(f"{n_warnings} instability warnings in log (>= {LINCS_WARN_COUNT_FOR_WARN})")

    energy_file = qc_dir / "energy.xvg"
    if energy_file.exists():
        data = read_xvg(energy_file)
        if data.size == 0:
            notes.append("energy.xvg empty")
            verdict = worst(verdict, "WARN")
        else:
            _, potential, total_e, temperature = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
            if not np.all(np.isfinite(potential)) or not np.all(np.isfinite(temperature)):
                verdict = worst(verdict, "FAIL")
                notes.append("non-finite (NaN/Inf) energy or temperature values")
            else:
                if temperature.max() > TEMP_MAX_K_FOR_FAIL:
                    verdict = worst(verdict, "FAIL")
                    notes.append(f"temperature spiked to {temperature.max():.0f} K (> {TEMP_MAX_K_FOR_FAIL:.0f} K)")
                diffs = np.abs(np.diff(potential))
                mad = np.median(np.abs(diffs - np.median(diffs))) or 1e-9
                spike_idx = np.where(diffs > ENERGY_SPIKE_MAD_MULTIPLE * mad)[0]
                if spike_idx.size:
                    t_spike = data[spike_idx[0] + 1, 0]
                    verdict = worst(verdict, "FAIL")
                    notes.append(f"potential-energy spike near t={t_spike:.1f} (frame-to-frame jump >> baseline noise)")
    else:
        notes.append("energy.xvg missing")
        verdict = worst(verdict, "WARN")

    return verdict, "; ".join(notes) if notes else "no warnings, energy/temperature stable"


def check_rmsd(qc_dir):
    f = qc_dir / "rmsd_backbone.xvg"
    if not f.exists():
        return "WARN", "rmsd_backbone.xvg missing"

    data = read_xvg(f)
    if data.shape[0] < 10:
        return "WARN", "too few RMSD frames to assess"

    t_ns, rmsd_nm = data[:, 0], data[:, 1]
    rmsd_a = rmsd_nm * 10.0

    n = len(rmsd_a)
    last_n = max(1, int(n * RMSD_LAST_FRAC))
    final_mean = rmsd_a[-last_n:].mean()

    half = n // 2
    slope_2nd, _ = np.polyfit(t_ns[half:], rmsd_a[half:], 1) if n - half >= 2 else (0.0, 0.0)

    verdict = "PASS"
    notes = [f"final {RMSD_LAST_FRAC*100:.0f}% mean RMSD = {final_mean:.2f} A",
             f"2nd-half slope = {slope_2nd:.4f} A/ns"]

    if final_mean > RMSD_FAIL_ANGSTROM or slope_2nd > RMSD_SECOND_HALF_SLOPE_FAIL:
        verdict = "FAIL"
    elif final_mean > RMSD_WARN_ANGSTROM or slope_2nd > RMSD_SECOND_HALF_SLOPE_WARN:
        verdict = "WARN"

    shape = "still drifting (not plateaued)" if slope_2nd > RMSD_SECOND_HALF_SLOPE_WARN else "plateaued / normal fluctuation"
    notes.append(shape)
    return verdict, "; ".join(notes)


def check_dissociation(qc_dir):
    mindist_f = qc_dir / "mindist_prot_lig.xvg"
    contacts_f = qc_dir / "numcontacts_prot_lig.xvg"
    if not mindist_f.exists() or not contacts_f.exists():
        return "WARN", "mindist/contacts files missing"

    mindist = read_xvg(mindist_f)
    contacts = read_xvg(contacts_f)
    if mindist.shape[0] < 2 or contacts.shape[0] < 2:
        return "WARN", "too few frames to assess"

    n = mindist.shape[0]
    last_n = max(1, int(n * DISSOC_LAST_FRAC))
    final_mindist_nm = mindist[-last_n:, 1].mean()
    final_ncontacts = contacts[-last_n:, 1].mean()

    notes = [f"final {DISSOC_LAST_FRAC*100:.0f}% mean min-dist = {final_mindist_nm*10:.2f} A "
             f"(contact cutoff {DISSOC_CONTACT_CUTOFF_NM*10:.1f} A), "
             f"mean contacts = {final_ncontacts:.1f}"]

    if final_mindist_nm > DISSOC_MINDIST_FAIL_NM and final_ncontacts < 0.5:
        return "FAIL", "peptide appears dissociated; " + notes[0]
    if final_mindist_nm > DISSOC_MINDIST_WARN_NM:
        return "WARN", "peptide loosely bound; " + notes[0]
    return "PASS", "peptide remains bound; " + notes[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=DEFAULT_ROOT)
    parser.add_argument("seeds", nargs="*", default=ALL_SEEDS)
    args = parser.parse_args()

    results = {}
    for seed in args.seeds:
        qc_dir = Path(args.root) / seed / "prod_md_48x1" / "qc"
        if not qc_dir.exists():
            results[seed] = {"clash": ("WARN", "qc/ directory missing -- run quality_check_extract.sh first"),
                              "rmsd": ("WARN", "n/a"), "dissoc": ("WARN", "n/a")}
            continue
        results[seed] = {
            "clash": check_clashes(qc_dir),
            "rmsd": check_rmsd(qc_dir),
            "dissoc": check_dissociation(qc_dir),
        }

    col = 18
    print(f"{'seed':<10} {'clashes/blowup':<{col}} {'RMSD trace':<{col}} {'peptide bound':<{col}} overall")
    print("-" * 100)
    any_fail = False
    for seed, r in results.items():
        overall = worst(r["clash"][0], r["rmsd"][0], r["dissoc"][0])
        any_fail = any_fail or overall == "FAIL"
        print(f"{seed:<10} {r['clash'][0]:<{col}} {r['rmsd'][0]:<{col}} {r['dissoc'][0]:<{col}} {overall}")

    print()
    for seed, r in results.items():
        print(f"--- {seed} ---")
        print(f"  clashes/blowup : [{r['clash'][0]}] {r['clash'][1]}")
        print(f"  RMSD trace     : [{r['rmsd'][0]}] {r['rmsd'][1]}")
        print(f"  peptide bound  : [{r['dissoc'][0]}] {r['dissoc'][1]}")

    sys.exit(1 if any_fail else 0)


if __name__ == "__main__":
    main()

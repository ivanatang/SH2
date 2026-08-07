#!/usr/bin/env python3
"""
Parse benchmark GROMACS logs from bchmk_prod_md_SH2.sh / submit_benchmarks_SH2.sh
and report ns/day performance per (mode, core count), averaged over replicates,
to identify the optimal core count and MPI strategy for this system.

Expected directory naming (under TARGET_ROOT):
  mpi_16cores_rep1/prod_bench.log
  thread_16cores_rep1/prod_bench.log
  ...

Usage:
  python3 parse_bchmk_performance_SH2.py [target_root]
  (defaults to /scratch/alpine/ivta1597/SH2/benchmark)
"""

from __future__ import annotations

import sys
from pathlib import Path
from collections import defaultdict

DEFAULT_TARGET_ROOT = Path("/scratch/alpine/ivta1597/SH2/benchmark")
LOG_NAME = "prod_bench.log"

import re

DIR_RE = re.compile(r"^(mpi|thread)_(\d+)cores_rep(\d+)$")
PERF_RE = re.compile(r"^\s*Performance:\s*([0-9]*\.?[0-9]+)\s+([0-9]*\.?[0-9]+)\s*$")


def extract_ns_per_day(log_path: Path) -> float:
    ns_day = None
    with log_path.open("r", errors="ignore") as f:
        for line in f:
            m = PERF_RE.match(line)
            if m:
                ns_day = float(m.group(1))
    if ns_day is None:
        raise ValueError(f"Could not find 'Performance:' line in {log_path}")
    return ns_day


def main():
    target_root = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_TARGET_ROOT
    if not target_root.exists():
        raise SystemExit(f"TARGET_ROOT does not exist: {target_root}")

    # (mode, cores) -> list of ns/day across reps
    results: dict[tuple[str, int], list[float]] = defaultdict(list)

    for p in sorted(target_root.iterdir()):
        if not p.is_dir():
            continue
        m = DIR_RE.match(p.name)
        if not m:
            continue
        mode, cores, rep = m.group(1), int(m.group(2)), int(m.group(3))

        log_path = p / LOG_NAME
        if not log_path.exists():
            print(f"[WARN] Missing log (job may still be running or failed): {log_path}")
            continue

        try:
            perf = extract_ns_per_day(log_path)
        except Exception as e:
            print(f"[WARN] Failed parsing {log_path}: {e}")
            continue

        results[(mode, cores)].append(perf)

    if not results:
        raise SystemExit(f"No matching benchmark directories/logs found under {target_root}")

    def mean(xs):
        return sum(xs) / len(xs)

    def std(xs):
        if len(xs) < 2:
            return 0.0
        m = mean(xs)
        return (sum((x - m) ** 2 for x in xs) / (len(xs) - 1)) ** 0.5

    rows = []
    for (mode, cores), perfs in results.items():
        rows.append((mode, cores, mean(perfs), std(perfs), len(perfs), perfs))

    # Sorted by performance, best first
    rows_by_perf = sorted(rows, key=lambda r: -r[2])

    print(f"TARGET_ROOT = {target_root}")
    print()
    print(f"{'mode':8s} {'cores':>6s} {'mean ns/day':>12s} {'std':>8s} {'n_reps':>7s}   raw values")
    print("-" * 80)
    for mode, cores, m, s, n, perfs in rows_by_perf:
        raw = ", ".join(f"{x:.2f}" for x in perfs)
        print(f"{mode:8s} {cores:6d} {m:12.2f} {s:8.2f} {n:7d}   [{raw}]")

    best_mode, best_cores, best_perf, *_ = rows_by_perf[0]
    print()
    print(f"BEST: mode={best_mode}, cores={best_cores}, {best_perf:.2f} ns/day (mean of {rows_by_perf[0][4]} rep(s))")

    # Also show, per mode, how performance scales with cores -- useful for
    # seeing where returns diminish (or go negative, like the original
    # catastrophic 64-rank real-MPI case) rather than just the single best point.
    print()
    print("Scaling by mode (sorted by core count):")
    for mode in sorted({r[0] for r in rows}):
        print(f"  {mode}:")
        mode_rows = sorted([r for r in rows if r[0] == mode], key=lambda r: r[1])
        for _, cores, m, s, n, _ in mode_rows:
            print(f"    {cores:3d} cores: {m:8.2f} +/- {s:5.2f} ns/day  (n={n})")


if __name__ == "__main__":
    main()

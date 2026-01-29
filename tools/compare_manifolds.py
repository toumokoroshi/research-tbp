#!/usr/bin/env python
"""
Compare two manifold output directories produced by StableManifoldComputation.

Usage:
  python tools/compare_manifolds.py <old_dir> <new_dir> [--max-files N] [--max-rows N]
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List, Tuple


def read_rows(path: Path, max_rows: int) -> List[List[float]]:
    rows: List[List[float]] = []
    with path.open("r", newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            try:
                vals = [float(x) for x in row]
            except ValueError:
                continue
            rows.append(vals)
            if max_rows and len(rows) >= max_rows:
                break
    return rows


def state_from_row(row: List[float]) -> Tuple[float, float, float, float, float, float]:
    # Expected: time,x,y,z,vx,vy,vz[,jacobi]
    return (row[1], row[2], row[3], row[4], row[5], row[6])


def jacobi_from_row(row: List[float]) -> float | None:
    return row[7] if len(row) > 7 else None


def diff_stats(a: List[List[float]], b: List[List[float]]) -> Dict[str, float]:
    n = min(len(a), len(b))
    if n == 0:
        return {
            "rows": 0,
            "max_pos": 0.0,
            "max_vel": 0.0,
            "max_state": 0.0,
            "mean_state": 0.0,
            "max_jacobi": 0.0,
        }

    max_pos = 0.0
    max_vel = 0.0
    max_state = 0.0
    sum_state = 0.0
    max_j = 0.0

    for i in range(n):
        sa = state_from_row(a[i])
        sb = state_from_row(b[i])
        dx = sa[0] - sb[0]
        dy = sa[1] - sb[1]
        dz = sa[2] - sb[2]
        dvx = sa[3] - sb[3]
        dvy = sa[4] - sb[4]
        dvz = sa[5] - sb[5]
        pos = math.sqrt(dx * dx + dy * dy + dz * dz)
        vel = math.sqrt(dvx * dvx + dvy * dvy + dvz * dvz)
        state = math.sqrt(pos * pos + vel * vel)
        max_pos = max(max_pos, pos)
        max_vel = max(max_vel, vel)
        max_state = max(max_state, state)
        sum_state += state

        ja = jacobi_from_row(a[i])
        jb = jacobi_from_row(b[i])
        if ja is not None and jb is not None:
            max_j = max(max_j, abs(ja - jb))

    return {
        "rows": float(n),
        "max_pos": max_pos,
        "max_vel": max_vel,
        "max_state": max_state,
        "mean_state": sum_state / n,
        "max_jacobi": max_j,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("old_dir", type=Path)
    parser.add_argument("new_dir", type=Path)
    parser.add_argument("--max-files", type=int, default=0, help="limit number of files")
    parser.add_argument("--max-rows", type=int, default=50, help="rows per file to compare (0=all)")
    args = parser.parse_args()

    old_dir = args.old_dir
    new_dir = args.new_dir
    if not old_dir.is_dir() or not new_dir.is_dir():
        raise SystemExit("Both old_dir and new_dir must be directories.")

    old_files = sorted(old_dir.glob("manifold_*.csv"))
    new_files = sorted(new_dir.glob("manifold_*.csv"))
    old_map = {p.name: p for p in old_files}
    new_map = {p.name: p for p in new_files}

    common = [name for name in old_map.keys() if name in new_map]
    common.sort()
    if args.max_files:
        common = common[: args.max_files]

    if not common:
        raise SystemExit("No common manifold_*.csv files found.")

    results = []
    for name in common:
        a = read_rows(old_map[name], args.max_rows)
        b = read_rows(new_map[name], args.max_rows)
        stats = diff_stats(a, b)
        results.append((name, stats))

    # Sort by max_state
    results.sort(key=lambda x: x[1]["max_state"], reverse=True)

    print(f"Compared {len(results)} files (max_rows={args.max_rows})")
    print("Top differences (max_state):")
    for name, s in results[:10]:
        print(
            f"  {name}: rows={int(s['rows'])} "
            f"max_state={s['max_state']:.3e} "
            f"max_pos={s['max_pos']:.3e} "
            f"max_vel={s['max_vel']:.3e} "
            f"mean_state={s['mean_state']:.3e} "
            f"max_jacobi={s['max_jacobi']:.3e}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

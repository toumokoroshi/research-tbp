#!/usr/bin/env python
"""
Plot Poincare section hits from manifold_section_hits.csv.

Usage:
  python tools/plot_manifold_section.py <csv_path> --out-dir <dir>
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt


def read_hits(path: Path) -> List[Dict[str, str]]:
    hits = []
    with path.open("r", newline="") as f:
        reader = csv.DictReader(
            filter(lambda row: row and not row.startswith("#"), f),
            fieldnames=[
                "manifold_idx",
                "type",
                "time",
                "x",
                "y",
                "z",
                "vx",
                "vy",
                "vz",
                "jacobi",
                "section_var",
                "section_value",
                "sign",
            ],
        )
        for row in reader:
            if not row or row["manifold_idx"].startswith("#"):
                continue
            hits.append(row)
    return hits


def to_float(hits: List[Dict[str, str]], key: str) -> List[float]:
    return [float(h[key]) for h in hits]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--out-dir", type=Path, default=Path("."))
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument("--split-type", action="store_true", help="separate stable/unstable plots")
    parser.add_argument(
        "--color-by-sign",
        action="store_true",
        help="color points by sign (direction) instead of type",
    )
    args = parser.parse_args()

    hits = read_hits(args.csv_path)
    if not hits:
        raise SystemExit("No hits found.")

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    components = ["x", "y", "z", "vx", "vy", "vz"]
    pairs = []
    for i in range(len(components)):
        for j in range(i + 1, len(components)):
            pairs.append((components[i], components[j]))

    types = sorted(set(h["type"] for h in hits))
    colors = {"stable": "tab:blue", "unstable": "tab:orange"}
    sign_colors = {"-1": "tab:purple", "1": "tab:green"}

    for a, b in pairs:
        if args.split_type:
            for t in types:
                subset = [h for h in hits if h["type"] == t]
                if not subset:
                    continue
                plt.figure(figsize=(6, 5))
                if args.color_by_sign:
                    for sign in ("-1", "1"):
                        sub2 = [h for h in subset if h["sign"] == sign]
                        xa = [float(h[a]) for h in sub2]
                        xb = [float(h[b]) for h in sub2]
                        if not xa:
                            continue
                        plt.scatter(
                            xa,
                            xb,
                            s=6,
                            alpha=0.6,
                            color=sign_colors.get(sign, "gray"),
                            label=f"sign {sign}",
                        )
                else:
                    xa = [float(h[a]) for h in subset]
                    xb = [float(h[b]) for h in subset]
                    plt.scatter(
                        xa,
                        xb,
                        s=6,
                        alpha=0.6,
                        color=colors.get(t, "gray"),
                        label=t,
                    )
                plt.xlabel(a)
                plt.ylabel(b)
                plt.title(f"Section Hits: {a} vs {b} ({t})")
                plt.legend()
                plt.tight_layout()
                out_path = out_dir / f"section_{a}_{b}_{t}.png"
                plt.savefig(out_path, dpi=args.dpi)
                plt.close()
        else:
            plt.figure(figsize=(6, 5))
            if args.color_by_sign:
                for sign in ("-1", "1"):
                    subset = [h for h in hits if h["sign"] == sign]
                    xa = [float(h[a]) for h in subset]
                    xb = [float(h[b]) for h in subset]
                    if not xa:
                        continue
                    plt.scatter(
                        xa,
                        xb,
                        s=6,
                        alpha=0.6,
                        color=sign_colors.get(sign, "gray"),
                        label=f"sign {sign}",
                    )
            else:
                for t in types:
                    subset = [h for h in hits if h["type"] == t]
                    xa = [float(h[a]) for h in subset]
                    xb = [float(h[b]) for h in subset]
                    plt.scatter(
                        xa,
                        xb,
                        s=6,
                        alpha=0.6,
                        color=colors.get(t, "gray"),
                        label=t,
                    )
            plt.xlabel(a)
            plt.ylabel(b)
            plt.title(f"Section Hits: {a} vs {b}")
            plt.legend()
            plt.tight_layout()
            out_path = out_dir / f"section_{a}_{b}.png"
            plt.savefig(out_path, dpi=args.dpi)
            plt.close()

    print(f"Saved {len(pairs)} plots to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

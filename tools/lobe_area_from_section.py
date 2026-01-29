#!/usr/bin/env python
"""
Compute lobe areas and plot lobe boundaries from manifold_section_hits.csv.
Alpha-shape is used by default for boundary extraction.

Usage:
  python tools/lobe_area_from_section.py <csv_path> --out-dir <dir>

Example:
  python tools/lobe_area_from_section.py data/.../manifold_section_hits.csv \
      --axes x vx --cj-bins 3.0007,3.0008,3.0009 --cluster-scale 4.0 --alpha-scale 1.5
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as mtri


@dataclass(frozen=True)
class Hit:
    manifold_idx: int
    type: str
    time: float
    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float
    jacobi: float
    section_var: str
    section_value: float
    sign: int

    def get(self, key: str) -> float:
        return getattr(self, key)


def read_hits(path: Path) -> List[Hit]:
    hits: List[Hit] = []
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
            hits.append(
                Hit(
                    manifold_idx=int(row["manifold_idx"]),
                    type=row["type"].strip(),
                    time=float(row["time"]),
                    x=float(row["x"]),
                    y=float(row["y"]),
                    z=float(row["z"]),
                    vx=float(row["vx"]),
                    vy=float(row["vy"]),
                    vz=float(row["vz"]),
                    jacobi=float(row["jacobi"]),
                    section_var=row["section_var"].strip(),
                    section_value=float(row["section_value"]),
                    sign=int(row["sign"]),
                )
            )
    return hits


def parse_axes(value: str) -> Tuple[str, str]:
    parts = value.split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("axes must be two comma-separated keys (e.g., x,vx)")
    return parts[0].strip(), parts[1].strip()


def parse_cj_bins(value: str) -> List[float]:
    parts = [p.strip() for p in value.split(",") if p.strip()]
    if len(parts) < 2:
        raise argparse.ArgumentTypeError("cj-bins must provide at least two edges")
    return [float(p) for p in parts]


def filter_hits(
    hits: Iterable[Hit],
    section_var: Optional[str],
    section_value: Optional[float],
    type_filter: Optional[str],
    sign_filter: Optional[int],
    cj_min: Optional[float],
    cj_max: Optional[float],
) -> List[Hit]:
    filtered = []
    for h in hits:
        if section_var and h.section_var != section_var:
            continue
        if section_value is not None and abs(h.section_value - section_value) > 1e-12:
            continue
        if type_filter and h.type != type_filter:
            continue
        if sign_filter is not None and h.sign != sign_filter:
            continue
        if cj_min is not None and h.jacobi < cj_min:
            continue
        if cj_max is not None and h.jacobi >= cj_max:
            continue
        filtered.append(h)
    return filtered


def estimate_spacing(points: Sequence[Tuple[float, float]]) -> float:
    if len(points) < 2:
        return 0.0
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    span_x = max(xs) - min(xs)
    span_y = max(ys) - min(ys)
    area = span_x * span_y
    if area <= 0.0:
        return 0.0
    return math.sqrt(area / len(points))


def build_grid(points: Sequence[Tuple[float, float]], cell_size: float) -> Dict[Tuple[int, int], List[int]]:
    grid: Dict[Tuple[int, int], List[int]] = {}
    inv = 1.0 / cell_size if cell_size > 0.0 else 1.0
    for i, (x, y) in enumerate(points):
        key = (int(math.floor(x * inv)), int(math.floor(y * inv)))
        grid.setdefault(key, []).append(i)
    return grid


def cluster_points(
    points: Sequence[Tuple[float, float]],
    eps: float,
) -> List[List[int]]:
    if len(points) == 0:
        return []
    if eps <= 0.0:
        return [list(range(len(points)))]

    grid = build_grid(points, eps)
    visited = [False] * len(points)
    clusters: List[List[int]] = []

    def neighbors(idx: int) -> List[int]:
        x, y = points[idx]
        key = (int(math.floor(x / eps)), int(math.floor(y / eps)))
        out: List[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                cell = (key[0] + dx, key[1] + dy)
                for j in grid.get(cell, []):
                    if j == idx:
                        continue
                    if math.hypot(points[j][0] - x, points[j][1] - y) <= eps:
                        out.append(j)
        return out

    for i in range(len(points)):
        if visited[i]:
            continue
        queue = [i]
        visited[i] = True
        cluster = [i]
        while queue:
            current = queue.pop()
            for j in neighbors(current):
                if not visited[j]:
                    visited[j] = True
                    queue.append(j)
                    cluster.append(j)
        clusters.append(cluster)
    return clusters


def order_points_angle(points: Sequence[Tuple[float, float]]) -> List[Tuple[float, float]]:
    if len(points) <= 2:
        return list(points)
    cx = sum(p[0] for p in points) / len(points)
    cy = sum(p[1] for p in points) / len(points)
    return sorted(points, key=lambda p: math.atan2(p[1] - cy, p[0] - cx))


def polygon_area(points: Sequence[Tuple[float, float]]) -> float:
    if len(points) < 3:
        return 0.0
    area = 0.0
    for i in range(len(points)):
        x1, y1 = points[i]
        x2, y2 = points[(i + 1) % len(points)]
        area += x1 * y2 - x2 * y1
    return abs(area) * 0.5


def build_lobes(
    hits: Sequence[Hit],
    axes: Tuple[str, str],
    cluster_scale: float,
    alpha_scale: float,
    use_alpha_shape: bool,
    min_points: int,
) -> List[Tuple[List[Tuple[float, float]], float]]:
    points = [(h.get(axes[0]), h.get(axes[1])) for h in hits]
    if not points:
        return []

    spacing = estimate_spacing(points)
    eps = spacing * cluster_scale if spacing > 0.0 else 0.0
    clusters = cluster_points(points, eps) if eps > 0.0 else [list(range(len(points)))]
    lobes: List[Tuple[List[Tuple[float, float]], float]] = []
    for cluster in clusters:
        if len(cluster) < min_points:
            continue
        cluster_points_list = [points[i] for i in cluster]
        if use_alpha_shape:
            lobes.extend(alpha_shape_lobes(cluster_points_list, alpha_scale, min_points))
        else:
            ordered = order_points_angle(cluster_points_list)
            area = polygon_area(ordered)
            lobes.append((ordered, area))
    return lobes


def alpha_shape_lobes(
    points: Sequence[Tuple[float, float]],
    alpha_scale: float,
    min_points: int,
) -> List[Tuple[List[Tuple[float, float]], float]]:
    if len(points) < 4:
        ordered = order_points_angle(points)
        area = polygon_area(ordered)
        return [(ordered, area)] if len(ordered) >= min_points else []

    spacing = estimate_spacing(points)
    alpha_radius = spacing * alpha_scale if spacing > 0.0 else 0.0
    if alpha_radius <= 0.0:
        ordered = order_points_angle(points)
        area = polygon_area(ordered)
        return [(ordered, area)] if len(ordered) >= min_points else []

    xs = np.array([p[0] for p in points], dtype=float)
    ys = np.array([p[1] for p in points], dtype=float)
    try:
        tri = mtri.Triangulation(xs, ys)
    except ValueError:
        ordered = order_points_angle(points)
        area = polygon_area(ordered)
        return [(ordered, area)] if len(ordered) >= min_points else []
    edges: Dict[Tuple[int, int], int] = {}
    for tri_indices in tri.triangles:
        i, j, k = int(tri_indices[0]), int(tri_indices[1]), int(tri_indices[2])
        ax, ay = xs[i], ys[i]
        bx, by = xs[j], ys[j]
        cx, cy = xs[k], ys[k]
        a = math.hypot(ax - bx, ay - by)
        b = math.hypot(bx - cx, by - cy)
        c = math.hypot(cx - ax, cy - ay)
        s = 0.5 * (a + b + c)
        area_sq = s * (s - a) * (s - b) * (s - c)
        if area_sq <= 0.0:
            continue
        circum_r = a * b * c / (4.0 * math.sqrt(area_sq))
        if circum_r > alpha_radius:
            continue
        for u, v in ((i, j), (j, k), (k, i)):
            key = (u, v) if u < v else (v, u)
            edges[key] = edges.get(key, 0) + 1

    boundary_edges = [e for e, count in edges.items() if count == 1]
    if not boundary_edges:
        return []

    adjacency: Dict[int, List[int]] = {}
    for u, v in boundary_edges:
        adjacency.setdefault(u, []).append(v)
        adjacency.setdefault(v, []).append(u)

    def edge_key(u: int, v: int) -> Tuple[int, int]:
        return (u, v) if u < v else (v, u)

    visited_edges: Dict[Tuple[int, int], bool] = {}
    lobes: List[Tuple[List[Tuple[float, float]], float]] = []
    for u, v in boundary_edges:
        if visited_edges.get(edge_key(u, v), False):
            continue
        path = [u, v]
        visited_edges[edge_key(u, v)] = True
        prev = u
        curr = v
        while True:
            neighbors = adjacency.get(curr, [])
            next_candidates = [n for n in neighbors if not visited_edges.get(edge_key(curr, n), False)]
            if not next_candidates:
                break
            nxt = next_candidates[0]
            visited_edges[edge_key(curr, nxt)] = True
            path.append(nxt)
            prev, curr = curr, nxt
            if curr == u:
                break

        if len(path) < min_points or path[-1] != path[0]:
            continue
        poly = [(float(xs[i]), float(ys[i])) for i in path]
        area = polygon_area(poly)
        if area > 0.0:
            lobes.append((poly, area))
    return lobes


def plot_lobes(
    out_path: Path,
    hits: Sequence[Hit],
    axes: Tuple[str, str],
    lobes_by_type: Dict[str, List[List[Tuple[float, float]]]],
    title: str,
    dpi: int,
) -> None:
    colors = {"stable": "tab:blue", "unstable": "tab:orange"}
    plt.figure(figsize=(7, 6))
    for t in sorted(set(h.type for h in hits)):
        subset = [h for h in hits if h.type == t]
        xa = [h.get(axes[0]) for h in subset]
        xb = [h.get(axes[1]) for h in subset]
        if xa:
            plt.scatter(xa, xb, s=6, alpha=0.6, color=colors.get(t, "gray"), label=t)
    for t, lobes in lobes_by_type.items():
        for poly in lobes:
            if len(poly) < 2:
                continue
            xs = [p[0] for p in poly] + [poly[0][0]]
            ys = [p[1] for p in poly] + [poly[0][1]]
            plt.plot(xs, ys, color=colors.get(t, "gray"), linewidth=1.2)
    plt.xlabel(axes[0])
    plt.ylabel(axes[1])
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi)
    plt.close()


def write_lobe_csv(
    out_path: Path,
    rows: List[Dict[str, object]],
) -> None:
    if not rows:
        return
    keys = list(rows[0].keys())
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--out-dir", type=Path, default=Path("."))
    parser.add_argument("--axes", type=parse_axes, default=("x", "vx"))
    parser.add_argument("--section-var", type=str, default=None)
    parser.add_argument("--section-value", type=float, default=None)
    parser.add_argument("--type", dest="type_filter", type=str, default=None)
    parser.add_argument("--sign", type=int, default=None)
    parser.add_argument("--cj-min", type=float, default=None)
    parser.add_argument("--cj-max", type=float, default=None)
    parser.add_argument("--cj-bins", type=parse_cj_bins, default=None)
    parser.add_argument("--cluster-scale", type=float, default=4.0)
    parser.add_argument("--alpha-scale", type=float, default=1.5)
    parser.add_argument("--use-alpha-shape", action="store_true", default=True)
    parser.add_argument("--no-alpha-shape", action="store_false", dest="use_alpha_shape")
    parser.add_argument("--min-points", type=int, default=10)
    parser.add_argument("--dpi", type=int, default=200)
    args = parser.parse_args()

    hits = read_hits(args.csv_path)
    if not hits:
        raise SystemExit("No hits found.")

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    axes = args.axes

    def process_band(cj_min: Optional[float], cj_max: Optional[float]) -> List[Dict[str, object]]:
        rows: List[Dict[str, object]] = []
        filtered = filter_hits(
            hits,
            section_var=args.section_var,
            section_value=args.section_value,
            type_filter=args.type_filter,
            sign_filter=args.sign,
            cj_min=cj_min,
            cj_max=cj_max,
        )
        if not filtered:
            return rows

        band_suffix = ""
        if cj_min is not None or cj_max is not None:
            band_suffix = f"_cj{cj_min:.6f}_{cj_max:.6f}" if cj_min is not None and cj_max is not None else "_cjband"
        for sign in sorted(set(h.sign for h in filtered)):
            sign_hits = [h for h in filtered if h.sign == sign]
            lobes_by_type: Dict[str, List[List[Tuple[float, float]]]] = {}
            for t in sorted(set(h.type for h in sign_hits)):
                subset = [h for h in sign_hits if h.type == t]
                lobes = build_lobes(
                    subset,
                    axes,
                    args.cluster_scale,
                    args.alpha_scale,
                    args.use_alpha_shape,
                    args.min_points,
                )
                lobes_by_type[t] = [poly for poly, _area in lobes]
                for idx, (poly, area) in enumerate(lobes):
                    rows.append(
                        {
                            "section_var": subset[0].section_var,
                            "section_value": subset[0].section_value,
                            "axis_a": axes[0],
                            "axis_b": axes[1],
                            "type": t,
                            "sign": sign,
                            "lobe_id": idx,
                            "area": area,
                            "cj_min": "" if cj_min is None else cj_min,
                            "cj_max": "" if cj_max is None else cj_max,
                            "num_points": len(poly),
                        }
                    )
            title = f"Lobe Boundaries: {axes[0]} vs {axes[1]} sign={sign}"
            plot_path = out_dir / f"lobe_section_{axes[0]}_{axes[1]}_sign{sign}{band_suffix}.png"
            plot_lobes(
                plot_path,
                sign_hits,
                axes,
                lobes_by_type,
                title,
                args.dpi,
            )
        return rows

    all_rows: List[Dict[str, object]] = []
    if args.cj_bins:
        for i in range(len(args.cj_bins) - 1):
            band_rows = process_band(args.cj_bins[i], args.cj_bins[i + 1])
            all_rows.extend(band_rows)
        out_path = out_dir / "lobe_areas_cj.csv"
    else:
        all_rows = process_band(args.cj_min, args.cj_max)
        out_path = out_dir / "lobe_areas.csv"

    if not all_rows:
        raise SystemExit("No lobes found. Try adjusting filters or cluster-scale.")

    write_lobe_csv(out_path, all_rows)
    print(f"Saved {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

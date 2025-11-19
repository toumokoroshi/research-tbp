#!/usr/bin/env python3
"""
Interactive SALI slice viewer.

Expected CSV header (case-insensitive):
Time,SALI,x,y,z,px,py,pz

- z is sliced with a slider.
- y is the vertical axis, x is the horizontal axis.
- Color shows the mean SALI in each x-y cell.
- Color map is anchored: SALI=-1→blue, SALI=0→white, SALI=1.41→red.
- Lines starting with '#' are treated as metadata comments.
- CSV can be chosen via --file-dialog at startup or the "Open CSV" button.

Extra quality-of-life options are available via CLI flags.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors
from matplotlib.widgets import Button, Slider


def _first_header_line(path: Path) -> str:
    """Return the first non-comment, non-empty line (lowercased) to detect the CSV header."""
    try:
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                return stripped.lower()
    except OSError:
        return ""
    return ""


def find_recent_sali_csv(search_root: Path) -> Optional[Path]:
    """Return the most recently modified CSV under search_root that looks like a SALI dump."""
    if not search_root.exists():
        return None

    newest: Optional[Path] = None
    for path in search_root.rglob("*.csv"):
        header = _first_header_line(path)
        if not header:
            continue
        if all(key in header for key in ("sali", "x", "y", "z")):
            if newest is None or path.stat().st_mtime > newest.stat().st_mtime:
                newest = path
    return newest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Slice SALI CSV data along z and plot x-y color contours with a z slider."
            " If no CSV is given, the newest SALI-like CSV under --search-root is used."
        )
    )
    parser.add_argument(
        "csv",
        nargs="?",
        type=Path,
        help="CSV file with columns Time,SALI,x,y,z,px,py,pz.",
    )
    parser.add_argument(
        "--file-dialog",
        action="store_true",
        help=(
            "Open a GUI file picker to choose the CSV (overrides the CLI path if selected)."
        ),
    )
    parser.add_argument(
        "--search-root",
        type=Path,
        default=Path("data"),
        help="Directory to search when CSV path is omitted (default: data).",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=120,
        help="Number of bins for the x-y grid (per axis).",
    )
    parser.add_argument(
        "--slice-width",
        type=float,
        default=0.02,
        help="Total thickness in z to include around the slider value.",
    )
    parser.add_argument(
        "--z",
        type=float,
        default=None,
        help="Initial z for the slice (default: median of z).",
    )
    parser.add_argument(
        "--xy-clip",
        type=float,
        default=0.0,
        help="Quantile for clipping x/y outliers (0 disables, typical: 0.01).",
    )
    parser.add_argument(
        "--sali-clip",
        type=float,
        default=0.001,
        help="Quantile for clipping SALI color scale (0 disables).",
    )
    parser.add_argument(
        "--downsample",
        type=int,
        default=0,
        help="Randomly keep at most this many rows (0 disables).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=120,
        help="Figure DPI.",
    )
    return parser.parse_args()


def load_data(path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    """Load CSV allowing for comment metadata lines that start with '#'."""
    metadata: dict[str, str] = {}
    try:
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith("#"):
                    stripped = stripped.lstrip("#").strip()
                    if "=" in stripped:
                        k, v = stripped.split("=", 1)
                        metadata[k.strip()] = v.strip()
                    continue
                # Stop reading once the header line is reached.
                break
    except OSError as exc:
        raise ValueError(f"Failed to read {path}: {exc}") from exc

    df = pd.read_csv(path, comment="#")
    df.columns = [c.strip().lower() for c in df.columns]

    expected = ("time", "sali", "x", "y", "z", "px", "py", "pz")
    missing = [col for col in expected if col not in df.columns]
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}")

    return df, metadata


def pick_file_via_dialog(initial_dir: Path) -> Optional[Path]:
    """Show a GUI file picker; return the chosen path or None if canceled/failure."""
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception:  # noqa: BLE001
        print("tkinter is unavailable; cannot open file dialog.", file=sys.stderr)
        return None

    root = tk.Tk()
    root.withdraw()
    # Avoid an empty/stale window on some platforms.
    root.update()
    filetypes = [("CSV files", "*.csv"), ("All files", "*.*")]
    path = filedialog.askopenfilename(
        title="Select SALI CSV",
        initialdir=str(initial_dir),
        filetypes=filetypes,
    )
    root.destroy()
    if not path:
        return None
    return Path(path)


def clip_xy(df: pd.DataFrame, quantile: float) -> pd.DataFrame:
    if quantile <= 0:
        return df
    quantile = min(max(quantile, 0.0), 0.3)
    lo = df[["x", "y"]].quantile(quantile)
    hi = df[["x", "y"]].quantile(1 - quantile)
    mask = df["x"].between(lo["x"], hi["x"]) & df["y"].between(lo["y"], hi["y"])
    return df.loc[mask]


def compute_limits(values: Iterable[float], clip_q: float) -> Tuple[float, float]:
    arr = np.asarray(list(values))
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return 0.0, 1.0
    if clip_q > 0:
        lo, hi = np.quantile(arr, [clip_q, 1 - clip_q])
    else:
        lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
    if lo == hi:
        hi = lo + 1e-12
    return float(lo), float(hi)


def make_mean_grid(
    df: pd.DataFrame,
    xedges: np.ndarray,
    yedges: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Return (x_centers, y_centers, mean Sali grid, sample count)."""
    x = df["x"].to_numpy()
    y = df["y"].to_numpy()
    w = df["sali"].to_numpy()

    hist_sum, _, _ = np.histogram2d(x, y, bins=[xedges, yedges], weights=w)
    hist_count, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])

    with np.errstate(invalid="ignore"):
        mean_grid = hist_sum / hist_count
    mean_grid = np.ma.array(mean_grid, mask=hist_count == 0)

    x_centers = 0.5 * (xedges[:-1] + xedges[1:])
    y_centers = 0.5 * (yedges[:-1] + yedges[1:])
    total_points = int(hist_count.sum())
    return x_centers, y_centers, mean_grid, total_points


def main() -> None:
    args = parse_args()

    csv_path = args.csv

    if args.file_dialog:
        chosen = pick_file_via_dialog(
            args.search_root if args.search_root else Path(".")
        )
        if chosen:
            csv_path = chosen
        else:
            print(
                "File dialog canceled; falling back to auto-detection.", file=sys.stderr
            )
    if csv_path is None:
        csv_path = find_recent_sali_csv(args.search_root)
        if csv_path is None:
            print(
                "CSV path not provided and no SALI-like CSV found under "
                f"{args.search_root}.",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"Auto-selected CSV: {csv_path}")

    if not csv_path.exists():
        print(f"CSV not found: {csv_path}", file=sys.stderr)
        sys.exit(1)

    try:
        df, metadata = load_data(csv_path)
    except Exception as exc:  # noqa: BLE001
        print(f"Failed to load CSV: {exc}", file=sys.stderr)
        sys.exit(1)

    if args.downsample and len(df) > args.downsample:
        df = df.sample(args.downsample, random_state=0).sort_index()
        print(f"Downsampled to {len(df):,} rows.")

    df = clip_xy(df, args.xy_clip).reset_index(drop=True)

    if args.z is None:
        z_init = float(df["z"].median())
    else:
        z_init = float(args.z)

    xedges = np.linspace(df["x"].min(), df["x"].max(), args.bins + 1)
    yedges = np.linspace(df["y"].min(), df["y"].max(), args.bins + 1)

    # Fixed anchors for SALI color mapping.
    sali_vmin = -1.0
    sali_vcenter = 0.0
    sali_vmax = 1.41
    sali_cmap = colors.LinearSegmentedColormap.from_list(
        "sali_anchor",
        [
            (0.0, "blue"),  # SALI=-1
            ((sali_vcenter - sali_vmin) / (sali_vmax - sali_vmin), "white"),  # SALI=0
            (1.0, "red"),  # SALI=1.41
        ],
    )

    fig, ax = plt.subplots(figsize=(9, 7), dpi=args.dpi)
    plt.subplots_adjust(left=0.22, bottom=0.24, right=0.88)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"SALI slice: {csv_path.name}")

    z_min = float(df["z"].min())
    z_max = float(df["z"].max())
    z_span = max(z_max - z_min, 1e-6)
    width_min = max(z_span / 2000, 1e-6)
    width_max = max(z_span, width_min)

    slider_ax = fig.add_axes([0.22, 0.13, 0.65, 0.03])
    width_ax = fig.add_axes([0.22, 0.08, 0.65, 0.03])
    save_ax = fig.add_axes([0.05, 0.36, 0.12, 0.05])
    file_ax = fig.add_axes([0.05, 0.28, 0.12, 0.05])

    z_slider = Slider(
        slider_ax,
        "z slice",
        valmin=z_min,
        valmax=z_max,
        valinit=z_init,
        valstep=z_span / 500,
    )
    width_slider = Slider(
        width_ax,
        "z width",
        valmin=width_min,
        valmax=width_max,
        valinit=min(max(args.slice_width, width_min), width_max),
    )

    file_button = Button(file_ax, "Open CSV", color="#ddeeff", hovercolor="#bbccff")
    save_button = Button(save_ax, "Save PNG", color="#dddddd", hovercolor="#bbbbbb")

    info_text = ax.text(
        0.02,
        0.98,
        "",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.65),
    )

    contour = None
    colorbar = None

    def meta_label_text(meta: dict[str, str]) -> str:
        if not meta:
            return ""
        parts = []
        if "MESH CENTER" in meta:
            parts.append(f"center: {meta['MESH CENTER']}")
        if "MESH DIVISION" in meta:
            parts.append(f"div: {meta['MESH DIVISION']}")
        if "MESH HALF WIDTH" in meta:
            parts.append(f"half width: {meta['MESH HALF WIDTH']}")
        if "SIMULATION TIME" in meta:
            parts.append(f"T={meta['SIMULATION TIME']}")
        if "CALCULATION TIMESTEP" in meta:
            parts.append(f"dt={meta['CALCULATION TIMESTEP']}")
        return "\n".join(parts)

    meta_text = ax.text(
        0.98,
        0.02,
        meta_label_text(metadata),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.55),
    )

    def set_data(new_path: Path) -> bool:
        nonlocal df, csv_path, xedges, yedges, metadata
        try:
            new_df, new_meta = load_data(new_path)
        except Exception as exc:  # noqa: BLE001
            print(f"Failed to load CSV: {exc}", file=sys.stderr)
            return False

        if args.downsample and len(new_df) > args.downsample:
            new_df = new_df.sample(args.downsample, random_state=0).sort_index()
            print(f"Downsampled to {len(new_df):,} rows.")

        new_df = clip_xy(new_df, args.xy_clip).reset_index(drop=True)
        if new_df.empty:
            print(f"{new_path} has no rows after clipping.", file=sys.stderr)
            return False

        csv_path = new_path
        df = new_df
        metadata = new_meta.copy()
        xedges = np.linspace(df["x"].min(), df["x"].max(), args.bins + 1)
        yedges = np.linspace(df["y"].min(), df["y"].max(), args.bins + 1)
        z_min_local = float(df["z"].min())
        z_max_local = float(df["z"].max())
        z_span_local = max(z_max_local - z_min_local, 1e-6)
        z_slider.valmin = z_min_local
        z_slider.valmax = z_max_local
        z_slider.valstep = z_span_local / 500
        z_slider.ax.set_xlim(z_min_local, z_max_local)
        z_slider.set_val(np.clip(float(df["z"].median()), z_min_local, z_max_local))

        width_min_local = max(z_span_local / 2000, 1e-6)
        width_max_local = max(z_span_local, width_min_local)
        width_slider.valmin = width_min_local
        width_slider.valmax = width_max_local
        width_slider.ax.set_xlim(width_min_local, width_max_local)
        width_slider.set_val(
            min(max(args.slice_width, width_min_local), width_max_local)
        )

        ax.set_title(f"SALI slice: {csv_path.name}")
        meta_text.set_text(meta_label_text(metadata))
        redraw_slice()
        return True

    def current_norm_and_levels() -> Tuple[colors.Normalize, np.ndarray]:
        norm = colors.TwoSlopeNorm(vmin=sali_vmin, vcenter=sali_vcenter, vmax=sali_vmax)
        levels = np.linspace(sali_vmin, sali_vmax, 48)
        return norm, levels

    def redraw_slice(_=None) -> None:
        nonlocal contour, colorbar
        z_center = z_slider.val
        half_width = width_slider.val * 0.5
        mask = (df["z"] >= z_center - half_width) & (df["z"] <= z_center + half_width)
        sliced = df.loc[mask]

        # Clear previous contour artists from the axes.
        for artist in list(ax.collections):
            try:
                artist.remove()
            except Exception:
                pass

        if sliced.empty:
            info_text.set_text(
                f"z={z_center:.4g} +/- {half_width:.4g}\n(no data in slice)"
            )
            ax.figure.canvas.draw_idle()
            return

        x_centers, y_centers, grid, sample_count = make_mean_grid(
            sliced, xedges, yedges
        )
        if grid.mask.all():
            info_text.set_text(
                f"z={z_center:.4g} +/- {half_width:.4g}\n(no data in grid)"
            )
            ax.figure.canvas.draw_idle()
            return

        Xc, Yc = np.meshgrid(x_centers, y_centers, indexing="ij")
        norm, levels = current_norm_and_levels()
        contour = ax.contourf(
            Xc,
            Yc,
            grid,
            levels=levels,
            cmap=sali_cmap,
            norm=norm,
            extend="both",
        )
        if colorbar is None:
            colorbar = fig.colorbar(contour, ax=ax)
            colorbar.set_label("mean SALI in cell")
        else:
            colorbar.update_normal(contour)

        info_text.set_text(
            f"z={z_center:.4g} +/- {half_width:.4g}\n"
            f"points={len(sliced):,} (grid uses {sample_count:,})"
        )

        ax.set_xlim(xedges[0], xedges[-1])
        ax.set_ylim(yedges[0], yedges[-1])
        ax.figure.canvas.draw_idle()

    def on_open(_event) -> None:
        start_dir = csv_path.parent if csv_path else args.search_root
        chosen = pick_file_via_dialog(start_dir)
        if chosen and chosen.exists():
            if set_data(chosen):
                print(f"Loaded {chosen}")
        elif chosen:
            print(f"Selected file does not exist: {chosen}", file=sys.stderr)

    def on_save(_event) -> None:
        z_center = z_slider.val
        half_width = width_slider.val * 0.5
        out_name = f"{csv_path.stem}_z{z_center:+.3f}_w{2*half_width:.3f}_sali.png"
        out_path = csv_path.parent / out_name
        fig.savefig(out_path, dpi=max(args.dpi, 150))
        print(f"Saved current view to {out_path}")

    z_slider.on_changed(redraw_slice)
    width_slider.on_changed(redraw_slice)
    file_button.on_clicked(on_open)
    save_button.on_clicked(on_save)

    redraw_slice()
    plt.show()


if __name__ == "__main__":
    main()

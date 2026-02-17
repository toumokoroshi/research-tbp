"""SALI時系列解析スクリプト.

sali_inertial_test の出力CSVを読み込み、
log(SALI) vs log(t) と log(SALI) vs t をプロットし、
べき乗則（規則軌道・周期軌道）vs 指数関数（カオス軌道）を判別する。

入力ファイル形式: time,sali,w1_w2_dot,qx,qy,qz (CSV)
"""

from pathlib import Path
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

try:
    import japanize_matplotlib  # noqa: F401
except ImportError:
    pass

# --- 設定 ---
FIGURE_DPI: int = 300


def load_sali_timeseries(filepath: Path) -> dict[str, np.ndarray]:
    """SALI時系列CSVを読み込む。

    Args:
        filepath: CSVファイルのパス

    Returns:
        time, sali, w1_w2_dot, qx, qy, qz をキーとするnumpy配列の辞書
    """
    data = np.loadtxt(filepath, delimiter=",", skiprows=1)
    return {
        "time": data[:, 0],
        "sali": data[:, 1],
        "w1_w2_dot": data[:, 2],
        "qx": data[:, 3],
        "qy": data[:, 4],
        "qz": data[:, 5],
    }


def fit_power_law(t: np.ndarray, sali: np.ndarray) -> tuple[float, float]:
    """log-log空間で線形フィットし、べき乗指数を推定する。

    Args:
        t: 時刻配列（t > 0 のみ使用）
        sali: SALI配列（sali > 0 のみ使用）

    Returns:
        (exponent, intercept): SALI ~ C * t^exponent
    """
    mask = (t > 0) & (sali > 0)
    log_t = np.log10(t[mask])
    log_s = np.log10(sali[mask])
    if len(log_t) < 2:
        return 0.0, 0.0
    coeffs = np.polyfit(log_t, log_s, 1)
    return coeffs[0], coeffs[1]


def fit_exponential(t: np.ndarray, sali: np.ndarray) -> tuple[float, float]:
    """semi-log空間で線形フィットし、指数減衰率を推定する。

    Args:
        t: 時刻配列
        sali: SALI配列（sali > 0 のみ使用）

    Returns:
        (rate, intercept): SALI ~ C * exp(rate * t)
    """
    mask = sali > 0
    t_f = t[mask]
    log_s = np.log10(sali[mask])
    if len(t_f) < 2:
        return 0.0, 0.0
    coeffs = np.polyfit(t_f, log_s, 1)
    return coeffs[0], coeffs[1]


def plot_sali_analysis(
    datasets: dict[str, dict[str, np.ndarray]],
    title: str,
    output_path: Optional[Path] = None,
) -> None:
    """SALI時系列を4パネルでプロットする。

    Args:
        datasets: {ラベル: load_sali_timeseriesの戻り値} の辞書
        title: 図のタイトル
        output_path: 保存先パス（Noneなら表示のみ）
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    colors = ["#2196F3", "#FF5722", "#4CAF50", "#9C27B0"]

    for i, (label, data) in enumerate(datasets.items()):
        t = data["time"]
        sali = data["sali"]
        color = colors[i % len(colors)]

        mask_pos = (t > 0) & (sali > 0)

        # (1) log-log: べき乗則チェック
        ax = axes[0, 0]
        if mask_pos.any():
            ax.plot(np.log10(t[mask_pos]), np.log10(sali[mask_pos]),
                    label=label, color=color, alpha=0.8, linewidth=0.8)
            exponent, intercept = fit_power_law(t, sali)
            # フィットライン
            t_fit = np.linspace(np.log10(t[mask_pos].min()),
                                np.log10(t[mask_pos].max()), 100)
            ax.plot(t_fit, exponent * t_fit + intercept,
                    "--", color=color, alpha=0.5,
                    label=f"  fit: slope={exponent:.2f}")
        ax.set_xlabel("log₁₀(t)")
        ax.set_ylabel("log₁₀(SALI)")
        ax.set_title("log-log (べき乗則: 直線 → 規則/周期)")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        # (2) semi-log: 指数減衰チェック
        ax = axes[0, 1]
        if mask_pos.any():
            ax.plot(t[mask_pos], np.log10(sali[mask_pos]),
                    label=label, color=color, alpha=0.8, linewidth=0.8)
            rate, intercept_e = fit_exponential(t, sali)
            t_fit_lin = np.linspace(t[mask_pos].min(),
                                     t[mask_pos].max(), 100)
            ax.plot(t_fit_lin, rate * t_fit_lin + intercept_e,
                    "--", color=color, alpha=0.5,
                    label=f"  fit: rate={rate:.4f}")
        ax.set_xlabel("t")
        ax.set_ylabel("log₁₀(SALI)")
        ax.set_title("semi-log (指数減衰: 直線 → カオス)")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        # (3) SALI vs t（線形スケール）
        ax = axes[1, 0]
        ax.plot(t, sali, label=label, color=color, alpha=0.8, linewidth=0.8)
        ax.set_xlabel("t")
        ax.set_ylabel("SALI")
        ax.set_title("SALI(t) 線形スケール")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        # (4) w1·w2 内積（偏差ベクトルの整列度）
        ax = axes[1, 1]
        w1w2 = data["w1_w2_dot"]
        ax.plot(t, np.abs(w1w2), label=label, color=color,
                alpha=0.8, linewidth=0.8)
        ax.set_xlabel("t")
        ax.set_ylabel("|w1·w2|")
        ax.set_title("偏差ベクトル整列度 (1.0 = 完全整列)")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=14, fontweight="bold")
    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
        print(f"保存: {output_path}")
    plt.show()


def main() -> None:
    """メイン処理."""
    import argparse

    parser = argparse.ArgumentParser(description="SALI時系列解析")
    parser.add_argument(
        "input_dir",
        type=str,
        help="sali_inertial_test出力ディレクトリのパス",
    )
    parser.add_argument(
        "--orbit",
        type=str,
        default="LEO",
        choices=["LEO", "GEO", "both"],
        help="解析対象の軌道 (デフォルト: LEO)",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)

    orbits = ["LEO", "GEO"] if args.orbit == "both" else [args.orbit]

    for orbit in orbits:
        prefix = orbit
        files = {
            f"{orbit} 二体問題": input_dir / f"{prefix}_twobody.csv",
            f"{orbit} CR3BP(シンプレクティック)": input_dir / f"{prefix}_crtbp_symplectic.csv",
            f"{orbit} CR3BP(RK4)": input_dir / f"{prefix}_crtbp_rk.csv",
        }

        datasets: dict[str, dict[str, np.ndarray]] = {}
        for label, fpath in files.items():
            if fpath.exists():
                print(f"読み込み: {fpath}")
                datasets[label] = load_sali_timeseries(fpath)
            else:
                print(f"未検出（スキップ）: {fpath}")

        if datasets:
            output_path = input_dir / f"{orbit}_sali_analysis.png"
            plot_sali_analysis(
                datasets,
                title=f"{orbit} SALI時系列解析 — べき乗 vs 指数減衰",
                output_path=output_path,
            )

            # フィット結果サマリー
            print(f"\n{'='*60}")
            print(f"  {orbit} フィット結果サマリー")
            print(f"{'='*60}")
            for label, data in datasets.items():
                t = data["time"]
                sali = data["sali"]
                exp, _ = fit_power_law(t, sali)
                rate, _ = fit_exponential(t, sali)
                print(f"  {label}:")
                print(f"    べき乗フィット: slope = {exp:.3f}"
                      f"  (周期軌道なら ≈ -2.0)")
                print(f"    指数フィット:   rate  = {rate:.6f}"
                      f"  (カオスなら大きな負値)")
                final_sali = sali[sali > 0][-1] if (sali > 0).any() else 0
                print(f"    最終SALI: {final_sali:.6e}")
            print()


if __name__ == "__main__":
    main()

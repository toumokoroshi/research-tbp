"""
Zero Velocity Curves (ZVC) and Forbidden Regions for the Sun-Earth System (2D)

Interactive visualization with adjustable Jacobi integral values and zoom.
Displays Lagrange points L1-L5.

Author: Research Team
Date: 2026-01-29
"""

from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button, TextBox

# Font settings
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"

# ZVC Style (unified with saliplotter.py)
ZERO_VELOCITY_COLOR = "#d35400"
FORBIDDEN_REGION_COLOR = "#ffe8d6"
# Sun-Earth system mass parameter μ = M_earth / (M_sun + M_earth)
MU_SUN_EARTH: float = 3.0034805945421924e-06
SLIDER_RES_X: float = 0.010000000000000
SLIDER_RES_Y: float = 0.010000000000000000
SLIDER_RES_ZOOM: float = 0.00100000000000000000000
SLIDER_RES_JACOBI: float = 0.00001000000000000000000000
SLIDER_RES_FIGSIZE: float = 1.000000000000000000000000000000
SLIDER_RES_DPI: int = 10


def compute_effective_potential(
    x: np.ndarray, y: np.ndarray, mu: float
) -> np.ndarray:
    """Compute the effective potential (Jacobi integral C = 2Ω)."""
    x1: float = -mu
    x2: float = 1.0 - mu

    r1: np.ndarray = np.sqrt((x - x1) ** 2 + y ** 2)
    r2: np.ndarray = np.sqrt((x - x2) ** 2 + y ** 2)

    r1 = np.maximum(r1, 1e-10)
    r2 = np.maximum(r2, 1e-10)

    omega: np.ndarray = 0.5 * (x ** 2 + y ** 2) + (1 - mu) / r1 + mu / r2
    return 2 * omega


def compute_lagrange_points(mu: float) -> dict[str, Tuple[float, float]]:
    """Compute Lagrange points L1-L5."""
    from scipy.optimize import brentq

    x2: float = 1.0 - mu

    def eq_l1(x: float) -> float:
        return x - (1 - mu) / (x + mu) ** 2 + mu / (x - 1 + mu) ** 2

    def eq_l2(x: float) -> float:
        return x - (1 - mu) / (x + mu) ** 2 - mu / (x - 1 + mu) ** 2

    def eq_l3(x: float) -> float:
        return x + (1 - mu) / (x + mu) ** 2 + mu / (x - 1 + mu) ** 2

    x_l1: float = brentq(eq_l1, 0.5, x2 - 0.001)
    x_l2: float = brentq(eq_l2, x2 + 0.001, 1.5)
    x_l3: float = brentq(eq_l3, -1.5, -mu - 0.001)

    x_l4: float = 0.5 - mu
    y_l4: float = np.sqrt(3) / 2
    y_l5: float = -np.sqrt(3) / 2

    return {
        "L1": (x_l1, 0.0),
        "L2": (x_l2, 0.0),
        "L3": (x_l3, 0.0),
        "L4": (x_l4, y_l4),
        "L5": (x_l4, y_l5),
    }


def compute_jacobi_at_lagrange_points(
    lagrange_points: dict[str, Tuple[float, float]], mu: float
) -> dict[str, float]:
    """Compute Jacobi integral values at each Lagrange point."""
    jacobi_values: dict[str, float] = {}
    for name, (x, y) in lagrange_points.items():
        c = compute_effective_potential(np.array([x]), np.array([y]), mu)
        jacobi_values[name] = float(c[0])
    return jacobi_values


class ZVCPlotter2D:
    """Interactive 2D Zero Velocity Curve and Forbidden Region plotter."""

    def __init__(
        self,
        mu: float = MU_SUN_EARTH,
        resolution: int = 3000,
        preview_resolution: int = 500,
    ) -> None:
        self.mu = mu
        self.resolution = resolution  # High-res for saving
        self.preview_resolution = preview_resolution  # Low-res for interactive

        self.lagrange_points = compute_lagrange_points(mu)
        self.jacobi_at_lp = compute_jacobi_at_lagrange_points(
            self.lagrange_points, mu
        )

        # Current values
        self.current_jacobi = self.jacobi_at_lp["L1"]
        self.current_xcenter = 1.0
        self.current_ycenter = 0.0
        self.current_zoom = 0.03

    def _compute_grid_for_view(
        self, xlim: Tuple[float, float], ylim: Tuple[float, float],
        use_high_res: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute grid dynamically for current view.
        
        Args:
            xlim: X軸の範囲
            ylim: Y軸の範囲
            use_high_res: Trueなら高解像度(保存用)、Falseなら低解像度(プレビュー用)
        """
        res = self.resolution if use_high_res else self.preview_resolution
        margin = 0.1 * max(xlim[1] - xlim[0], ylim[1] - ylim[0])
        x = np.linspace(xlim[0] - margin, xlim[1] + margin, res)
        y = np.linspace(ylim[0] - margin, ylim[1] + margin, res)
        X, Y = np.meshgrid(x, y)
        C = compute_effective_potential(X, Y, self.mu)
        return X, Y, C

    def plot(self, ax: plt.Axes, use_high_res: bool = False) -> None:
        """Plot ZVC and forbidden regions.
        
        Args:
            ax: matplotlib Axes
            use_high_res: Trueなら高解像度(保存用)、Falseなら低解像度(プレビュー用)
        """
        ax.clear()

        half_size = self.current_zoom / 2
        xlim = (self.current_xcenter - half_size, self.current_xcenter + half_size)
        ylim = (self.current_ycenter - half_size, self.current_ycenter + half_size)

        X, Y, C = self._compute_grid_for_view(xlim, ylim, use_high_res=use_high_res)

        # Forbidden region: where C < jacobi_c
        c_min = C.min()
        if c_min < self.current_jacobi:
            ax.contourf(
                X, Y, C,
                levels=[c_min, self.current_jacobi],
                colors=[FORBIDDEN_REGION_COLOR],
                alpha=0.7,
            )

        # Zero velocity curve (only if jacobi_c is within the data range)
        if c_min <= self.current_jacobi <= C.max():
            ax.contour(
                X, Y, C,
                levels=[self.current_jacobi],
                colors=[ZERO_VELOCITY_COLOR],
                linewidths=2.0,
            )

        # Sun and Earth positions
        x_sun = -self.mu
        x_earth = 1.0 - self.mu
        # ax.plot(x_sun, 0, "yo", markersize=15, label="Sun")
        # ax.plot(x_earth, 0, "bo", markersize=8, label="Earth")
        ax.plot(x_sun, 0, "yo", markersize=15)
        ax.plot(x_earth, 0, "bo", markersize=8)

        # Earth annotation
        ax.annotate(
            "EARTH", (x_earth, 0),
            textcoords="offset points",
            xytext=(0, -20),
            fontsize=12,
            fontweight="bold",
            color="darkblue",
            ha="center",
        )
        # sun annotation
        ax.annotate(
            "SUN", (x_sun, 0),
            textcoords="offset points",
            xytext=(0, -20),
            fontsize=12,
            fontweight="bold",
            color="darkorange",
            ha="center",
        )

        # Forbidden region annotations
        # 禁止領域の上部と下部の中央にラベルを配置
        half_size = self.current_zoom / 2
        label_y_top = self.current_ycenter + half_size * 0.7
        label_y_bottom = self.current_ycenter - half_size * 0.7
        label_x = self.current_xcenter

        # ax.text(
        #     label_x, label_y_top, "FORBIDDEN REGION",
        #     fontsize=14, fontweight="bold", color=ZERO_VELOCITY_COLOR,
        #     ha="center", va="center",
        # )
        # ax.text(
        #     label_x, label_y_bottom, "FORBIDDEN REGION",
        #     fontsize=14, fontweight="bold", color=ZERO_VELOCITY_COLOR,
        #     ha="center", va="center",
        # )

        # Lagrange points (with arrow lines instead of markers)
        for name, (x, y) in self.lagrange_points.items():
            # 矢印の方向とオフセットを設定
            if name in ["L1"]:
                xytext = (-30, 15)
            elif name in ["L2", "L3"]:
                xytext = (30, 15)
            elif name == "L4":
                # L4は右上方向
                xytext = (20, 15)
            else:  # L5
                # L5は右下方向
                xytext = (20, -15)
            
            ax.annotate(
                name, (x, y),
                textcoords="offset points",
                xytext=xytext,
                fontsize=10,
                color="#8e44ad",
                ha="center",
                arrowprops=dict(
                    arrowstyle="-",
                    color="#8e44ad",
                    lw=1.0,
                ),
            )

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal")
        ax.set_xlabel("x [AU]", fontsize=12)
        ax.set_ylabel("y [AU]", fontsize=12)
        
        # 軸のオフセット表記を無効化（絶対値で表示）
        ax.ticklabel_format(useOffset=False, style='plain')
        
        ax.set_title(
            f"ZVC & Forbidden Regions\n"
            f"Jacobi C = {self.current_jacobi:.10f}",
            fontsize=14,
        )
        # ax.legend(loc="upper right")
        ax.grid(True, linestyle="--", alpha=0.5)

        # Lagrange point info
        info_text = "C at L-points:\n"
        for name in ["L1", "L2", "L3", "L4", "L5"]:
            info_text += f"{name}: {self.jacobi_at_lp[name]:.10f}\n"
        # ax.text(
        #     0.02, 0.98, info_text,
        #     transform=ax.transAxes,
        #     fontsize=9,
        #     verticalalignment="top",
        #     bbox={"boxstyle": "round", "facecolor": "wheat", "alpha": 0.8},
        # )

    def run_interactive(self) -> None:
        """Run in interactive mode with zoom and pan controls."""
        fig, ax = plt.subplots(figsize=(10, 10))
        plt.subplots_adjust(bottom=0.50, right=0.85)

        # Image save settings
        self.save_figsize = 10
        self.save_dpi = 300

        # Sliders
        ax_jacobi = plt.axes([0.15, 0.30, 0.50, 0.03])
        slider_jacobi = Slider(
            ax_jacobi, "Jacobi C", 2.99, 3.05,
            valinit=self.current_jacobi, valstep=SLIDER_RES_JACOBI,
        )

        ax_xcenter = plt.axes([0.15, 0.24, 0.50, 0.03])
        slider_xcenter = Slider(
            ax_xcenter, "X Center", -2.0, 2.0,
            valinit=0.0, valstep=SLIDER_RES_X,
        )

        ax_ycenter = plt.axes([0.15, 0.18, 0.50, 0.03])
        slider_ycenter = Slider(
            ax_ycenter, "Y Center", -2.0, 2.0, 
            valinit=0.0, valstep=SLIDER_RES_Y,
        )

        ax_zoom = plt.axes([0.15, 0.12, 0.50, 0.03])
        slider_zoom = Slider(
            ax_zoom, "View Size", 0.0001, 4.0,
            valinit=3.0, valstep=SLIDER_RES_ZOOM,
        )

        # TextBoxes for keyboard input
        ax_tb_jacobi = plt.axes([0.75, 0.30, 0.12, 0.03])
        tb_jacobi = TextBox(ax_tb_jacobi, "", initial=f"{self.current_jacobi:.10f}")

        ax_tb_xcenter = plt.axes([0.75, 0.24, 0.12, 0.03])
        tb_xcenter = TextBox(ax_tb_xcenter, "", initial="1.0")

        ax_tb_ycenter = plt.axes([0.75, 0.18, 0.12, 0.03])
        tb_ycenter = TextBox(ax_tb_ycenter, "", initial="0.0")

        ax_tb_zoom = plt.axes([0.75, 0.12, 0.12, 0.03])
        tb_zoom = TextBox(ax_tb_zoom, "", initial="0.03")

        # Image save size settings
        ax_figsize = plt.axes([0.15, 0.06, 0.50, 0.03])
        slider_figsize = Slider(
            ax_figsize, "Fig Size", 4, 20,
            valinit=self.save_figsize, valstep=SLIDER_RES_FIGSIZE,
        )
        ax_tb_figsize = plt.axes([0.75, 0.06, 0.12, 0.03])
        tb_figsize = TextBox(ax_tb_figsize, "", initial=f"{self.save_figsize}")

        ax_dpi = plt.axes([0.15, 0.01, 0.50, 0.03])
        slider_dpi = Slider(
            ax_dpi, "DPI", 72, 600,
            valinit=self.save_dpi, valstep=SLIDER_RES_DPI,
        )
        ax_tb_dpi = plt.axes([0.75, 0.01, 0.12, 0.03])
        tb_dpi = TextBox(ax_tb_dpi, "", initial=f"{self.save_dpi}")

        # Size info label
        ax_info = plt.axes([0.15, 0.42, 0.70, 0.03])
        ax_info.axis('off')
        size_info = ax_info.text(
            0.0, 0.5, 
            f"Output: {self.save_figsize * self.save_dpi} x {self.save_figsize * self.save_dpi} px",
            fontsize=10, va='center'
        )

        # Save button
        ax_save = plt.axes([0.88, 0.01, 0.10, 0.05])
        btn_save = Button(ax_save, "Save")

        def update_from_sliders(_) -> None:
            self.current_jacobi = slider_jacobi.val
            self.current_xcenter = slider_xcenter.val
            self.current_ycenter = slider_ycenter.val
            self.current_zoom = slider_zoom.val
            # Update textboxes
            tb_jacobi.set_val(f"{self.current_jacobi:.10f}")
            tb_xcenter.set_val(f"{self.current_xcenter:.4f}")
            tb_ycenter.set_val(f"{self.current_ycenter:.4f}")
            tb_zoom.set_val(f"{self.current_zoom:.4f}")
            self.plot(ax)
            fig.canvas.draw_idle()

        def update_jacobi_from_text(text: str) -> None:
            try:
                val = float(text)
                if 2.99 <= val <= 3.05:
                    self.current_jacobi = val
                    slider_jacobi.set_val(val)
            except ValueError:
                pass

        def update_xcenter_from_text(text: str) -> None:
            try:
                val = float(text)
                self.current_xcenter = val
                slider_xcenter.set_val(np.clip(val, -2.0, 2.0))
                self.plot(ax)
                fig.canvas.draw_idle()
            except ValueError:
                pass

        def update_ycenter_from_text(text: str) -> None:
            try:
                val = float(text)
                self.current_ycenter = val
                slider_ycenter.set_val(np.clip(val, -2.0, 2.0))
                self.plot(ax)
                fig.canvas.draw_idle()
            except ValueError:
                pass

        def update_zoom_from_text(text: str) -> None:
            try:
                val = float(text)
                if val > 0:
                    self.current_zoom = val
                    slider_zoom.set_val(np.clip(val, 0.001, 4.0))
                    self.plot(ax)
                    fig.canvas.draw_idle()
            except ValueError:
                pass

        def update_figsize_from_slider(_) -> None:
            self.save_figsize = int(slider_figsize.val)
            tb_figsize.set_val(f"{self.save_figsize}")
            update_size_info()

        def update_dpi_from_slider(_) -> None:
            self.save_dpi = int(slider_dpi.val)
            tb_dpi.set_val(f"{self.save_dpi}")
            update_size_info()

        def update_figsize_from_text(text: str) -> None:
            try:
                val = int(float(text))
                if 4 <= val <= 20:
                    self.save_figsize = val
                    slider_figsize.set_val(val)
                    update_size_info()
            except ValueError:
                pass

        def update_dpi_from_text(text: str) -> None:
            try:
                val = int(float(text))
                if 72 <= val <= 600:
                    self.save_dpi = val
                    slider_dpi.set_val(val)
                    update_size_info()
            except ValueError:
                pass

        def update_size_info() -> None:
            px = self.save_figsize * self.save_dpi
            size_info.set_text(f"Output: {px} x {px} px")
            fig.canvas.draw_idle()

        def save_image(_) -> None:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            default_name = f"zvc_C{self.current_jacobi:.10f}_x{self.current_xcenter:.2f}_y{self.current_ycenter:.2f}_z{self.current_zoom:.2f}.jpg"
            filepath = filedialog.asksaveasfilename(
                defaultextension=".jpg",
                filetypes=[("JPEG", "*.jpg"), ("PNG", "*.png"), ("All files", "*.*")],
                initialfile=default_name,
                title="Save ZVC Image"
            )
            root.destroy()
            if filepath:
                save_fig, save_ax = plt.subplots(figsize=(self.save_figsize, self.save_figsize))
                self.plot(save_ax, use_high_res=True)  # High resolution for saving
                save_fig.savefig(filepath, dpi=self.save_dpi, bbox_inches="tight")
                plt.close(save_fig)
                px = self.save_figsize * self.save_dpi
                print(f"Saved: {filepath} ({px}x{px}px)")

        slider_jacobi.on_changed(update_from_sliders)
        slider_xcenter.on_changed(update_from_sliders)
        slider_ycenter.on_changed(update_from_sliders)
        slider_zoom.on_changed(update_from_sliders)

        tb_jacobi.on_submit(update_jacobi_from_text)
        tb_xcenter.on_submit(update_xcenter_from_text)
        tb_ycenter.on_submit(update_ycenter_from_text)
        tb_zoom.on_submit(update_zoom_from_text)

        slider_figsize.on_changed(update_figsize_from_slider)
        slider_dpi.on_changed(update_dpi_from_slider)
        tb_figsize.on_submit(update_figsize_from_text)
        tb_dpi.on_submit(update_dpi_from_text)

        btn_save.on_clicked(save_image)

        # Initial plot
        self.plot(ax)
        plt.show()

    def save_image(
        self,
        jacobi_c: float,
        x_center: float = 0.0,
        y_center: float = 0.0,
        zoom: float = 3.0,
        output_dir: Path = Path("."),
        filename: str | None = None,
    ) -> Path:
        """Save image for specified settings."""
        self.current_jacobi = jacobi_c
        self.current_xcenter = x_center
        self.current_ycenter = y_center
        self.current_zoom = zoom

        fig, ax = plt.subplots(figsize=(10, 10))
        self.plot(ax)
        if filename is None:
            filename = f"zvc_2d_C{jacobi_c:.10f}.jpg"
        save_path = output_dir / filename
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Image saved: {save_path}")
        return save_path


def main() -> None:
    """Main function."""
    plotter = ZVCPlotter2D()
    plotter.run_interactive()


if __name__ == "__main__":
    main()

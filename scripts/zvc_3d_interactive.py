"""
Zero Velocity Curves (ZVC) and Forbidden Regions for the Sun-Earth System (3D)

Interactive visualization with adjustable Jacobi integral values and zoom.
Displays 2D view and 3D ZVC isosurface with Lagrange points L1-L5.

Author: Research Team
Date: 2026-01-29
"""

from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button, TextBox
from mpl_toolkits.mplot3d import Axes3D
from skimage import measure

# Font settings
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"

# Sun-Earth system mass parameter μ = M_earth / (M_sun + M_earth)
MU_SUN_EARTH: float = 3.003467e-6
SLIDER_RES_X: float = 0.01
SLIDER_RES_Y: float = 0.01
SLIDER_RES_ZOOM: float = 0.001
SLIDER_RES_JACOBI: float = 0.00001


def compute_effective_potential_3d(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, mu: float
) -> np.ndarray:
    """Compute the 3D effective potential."""
    x1 = -mu
    x2 = 1.0 - mu

    r1 = np.sqrt((x - x1) ** 2 + y ** 2 + z ** 2)
    r2 = np.sqrt((x - x2) ** 2 + y ** 2 + z ** 2)

    r1 = np.maximum(r1, 1e-10)
    r2 = np.maximum(r2, 1e-10)

    omega = 0.5 * (x ** 2 + y ** 2) + (1 - mu) / r1 + mu / r2
    return 2 * omega


def compute_effective_potential_2d(
    x: np.ndarray, y: np.ndarray, mu: float
) -> np.ndarray:
    """Compute 2D effective potential (z=0 plane)."""
    return compute_effective_potential_3d(x, y, np.zeros_like(x), mu)


def compute_lagrange_points(mu: float) -> dict[str, Tuple[float, float, float]]:
    """Compute Lagrange points L1-L5."""
    from scipy.optimize import brentq

    x2 = 1.0 - mu

    def eq_l1(x):
        return x - (1 - mu) / (x + mu) ** 2 + mu / (x - 1 + mu) ** 2

    def eq_l2(x):
        return x - (1 - mu) / (x + mu) ** 2 - mu / (x - 1 + mu) ** 2

    def eq_l3(x):
        return x + (1 - mu) / (x + mu) ** 2 + mu / (x - 1 + mu) ** 2

    x_l1 = brentq(eq_l1, 0.5, x2 - 0.001)
    x_l2 = brentq(eq_l2, x2 + 0.001, 1.5)
    x_l3 = brentq(eq_l3, -1.5, -mu - 0.001)

    x_l4 = 0.5 - mu
    y_l4 = np.sqrt(3) / 2
    y_l5 = -np.sqrt(3) / 2

    return {
        "L1": (x_l1, 0.0, 0.0),
        "L2": (x_l2, 0.0, 0.0),
        "L3": (x_l3, 0.0, 0.0),
        "L4": (x_l4, y_l4, 0.0),
        "L5": (x_l4, y_l5, 0.0),
    }


def compute_jacobi_at_lagrange_points(
    lagrange_points: dict[str, Tuple[float, float, float]], mu: float
) -> dict[str, float]:
    """Compute Jacobi integral values at each Lagrange point."""
    jacobi_values = {}
    for name, (x, y, z) in lagrange_points.items():
        c = compute_effective_potential_3d(
            np.array([x]), np.array([y]), np.array([z]), mu
        )
        jacobi_values[name] = float(c[0])
    return jacobi_values


class ZVCPlotter3D:
    """3D Zero Velocity Surface plotter."""

    def __init__(self, mu: float = MU_SUN_EARTH, resolution_2d: int = 500, resolution_3d: int = 80) -> None:
        self.mu = mu
        self.resolution_2d = resolution_2d
        self.resolution_3d = resolution_3d

        self.lagrange_points = compute_lagrange_points(mu)
        self.jacobi_at_lp = compute_jacobi_at_lagrange_points(
            self.lagrange_points, mu
        )

        # Current values
        self.current_jacobi = self.jacobi_at_lp["L1"]
        self.current_xcenter = 0.0
        self.current_ycenter = 0.0
        self.current_zoom = 3.0

    def _compute_grid_2d(
        self, xlim: Tuple[float, float], ylim: Tuple[float, float]
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute 2D grid for current view."""
        margin = 0.1 * max(xlim[1] - xlim[0], ylim[1] - ylim[0])
        x = np.linspace(xlim[0] - margin, xlim[1] + margin, self.resolution_2d)
        y = np.linspace(ylim[0] - margin, ylim[1] + margin, self.resolution_2d)
        X, Y = np.meshgrid(x, y)
        C = compute_effective_potential_2d(X, Y, self.mu)
        return X, Y, C

    def _compute_zvc_isosurface(
        self, xlim: Tuple[float, float], ylim: Tuple[float, float], zlim: Tuple[float, float]
    ) -> Tuple[np.ndarray, np.ndarray] | Tuple[None, None]:
        """Compute ZVC isosurface using marching cubes."""
        margin = 0.05 * max(xlim[1] - xlim[0], ylim[1] - ylim[0])
        
        x = np.linspace(xlim[0] - margin, xlim[1] + margin, self.resolution_3d)
        y = np.linspace(ylim[0] - margin, ylim[1] + margin, self.resolution_3d)
        z = np.linspace(zlim[0], zlim[1], self.resolution_3d)
        
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        C = compute_effective_potential_3d(X, Y, Z, self.mu)
        
        # Use marching cubes to find the isosurface where C = jacobi_c
        try:
            verts, faces, _, _ = measure.marching_cubes(
                C, level=self.current_jacobi, 
                spacing=(x[1] - x[0], y[1] - y[0], z[1] - z[0])
            )
            # Shift vertices to actual coordinates
            verts[:, 0] += xlim[0] - margin
            verts[:, 1] += ylim[0] - margin
            verts[:, 2] += zlim[0]
            return verts, faces
        except ValueError:
            # No isosurface found
            return None, None

    def plot_2d(self, ax: plt.Axes) -> None:
        """Plot 2D ZVC and forbidden regions."""
        ax.clear()

        half = self.current_zoom / 2
        xlim = (self.current_xcenter - half, self.current_xcenter + half)
        ylim = (self.current_ycenter - half, self.current_ycenter + half)

        X, Y, C = self._compute_grid_2d(xlim, ylim)

        # Forbidden region: C < jacobi_c
        c_min = C.min()
        if c_min < self.current_jacobi:
            ax.contourf(
                X, Y, C,
                levels=[c_min, self.current_jacobi],
                colors=["gray"],
                alpha=1.0,
            )

        # Zero velocity curve
        if c_min <= self.current_jacobi <= C.max():
            ax.contour(
                X, Y, C,
                levels=[self.current_jacobi],
                colors="red",
                linewidths=2,
            )

        # Sun and Earth
        x_sun = -self.mu
        x_earth = 1.0 - self.mu
        ax.plot(x_sun, 0, "yo", markersize=12, label="Sun")
        ax.plot(x_earth, 0, "bo", markersize=6, label="Earth")

        # Lagrange points
        for name, (x, y, z) in self.lagrange_points.items():
            ax.plot(x, y, "r^", markersize=6)
            ax.annotate(
                name, (x, y),
                textcoords="offset points",
                xytext=(5, 5),
                fontsize=8,
                color="red",
            )

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal")
        ax.set_xlabel("x [AU]", fontsize=10)
        ax.set_ylabel("y [AU]", fontsize=10)
        ax.set_title(f"2D ZVC (z=0)\nC = {self.current_jacobi:.10f}", fontsize=11)
        ax.legend(loc="upper right", fontsize=7)
        ax.grid(True, linestyle="--", alpha=0.5)

    def plot_3d(self, ax3d) -> None:
        """Plot 3D ZVC isosurface."""
        ax3d.clear()

        half = self.current_zoom / 2
        xlim = (self.current_xcenter - half, self.current_xcenter + half)
        ylim = (self.current_ycenter - half, self.current_ycenter + half)
        zlim = (-half * 0.5, half * 0.5)  # z range is smaller

        # Compute ZVC isosurface
        verts, faces = self._compute_zvc_isosurface(xlim, ylim, zlim)
        
        if verts is not None and len(verts) > 0:
            # Plot the isosurface
            ax3d.plot_trisurf(
                verts[:, 0], verts[:, 1], verts[:, 2],
                triangles=faces,
                color="red",
                alpha=0.7,
                edgecolor="none",
            )

        # Sun and Earth on z=0 plane
        x_sun = -self.mu
        x_earth = 1.0 - self.mu
        ax3d.scatter([x_sun], [0], [0], color="yellow", s=100, label="Sun")
        ax3d.scatter([x_earth], [0], [0], color="blue", s=50, label="Earth")

        # Lagrange points
        for name, (lx, ly, lz) in self.lagrange_points.items():
            if xlim[0] <= lx <= xlim[1] and ylim[0] <= ly <= ylim[1]:
                ax3d.scatter([lx], [ly], [lz], color="green", s=40, marker="^")
                ax3d.text(lx, ly, lz + 0.05, name, fontsize=8, color="green")

        ax3d.set_xlim(xlim)
        ax3d.set_ylim(ylim)
        ax3d.set_zlim(zlim)
        ax3d.set_xlabel("x [AU]", fontsize=9)
        ax3d.set_ylabel("y [AU]", fontsize=9)
        ax3d.set_zlabel("z [AU]", fontsize=9)
        ax3d.set_title(f"3D ZVC Isosurface\nC = {self.current_jacobi:.10f}", fontsize=11)
        ax3d.legend(loc="upper right", fontsize=7)

    def run_interactive(self) -> None:
        """Run in interactive mode with zoom and pan controls."""
        fig = plt.figure(figsize=(16, 9))
        plt.subplots_adjust(bottom=0.35, left=0.05, right=0.95, wspace=0.25)

        # Create axes for plots
        ax2d = fig.add_axes([0.05, 0.40, 0.42, 0.55])
        ax3d = fig.add_axes([0.52, 0.40, 0.45, 0.55], projection="3d")

        # Sliders
        ax_jacobi = plt.axes([0.15, 0.25, 0.45, 0.03])
        slider_jacobi = Slider(
            ax_jacobi, "Jacobi C", 2.99, 3.05,
            valinit=self.current_jacobi, valstep=SLIDER_RES_JACOBI,
        )

        ax_xcenter = plt.axes([0.15, 0.19, 0.45, 0.03])
        slider_xcenter = Slider(
            ax_xcenter, "X Center", -2.0, 2.0,
            valinit=0.0, valstep=SLIDER_RES_X,
        )

        ax_ycenter = plt.axes([0.15, 0.13, 0.45, 0.03])
        slider_ycenter = Slider(
            ax_ycenter, "Y Center", -2.0, 2.0,
            valinit=0.0, valstep=SLIDER_RES_Y,
        )

        ax_zoom = plt.axes([0.15, 0.07, 0.45, 0.03])
        slider_zoom = Slider(
            ax_zoom, "View Size", 0.001, 4.0,
            valinit=3.0, valstep=SLIDER_RES_ZOOM,
        )

        # TextBoxes for keyboard input
        ax_tb_jacobi = plt.axes([0.70, 0.25, 0.12, 0.03])
        tb_jacobi = TextBox(ax_tb_jacobi, "", initial=f"{self.current_jacobi:.10f}")

        ax_tb_xcenter = plt.axes([0.70, 0.19, 0.12, 0.03])
        tb_xcenter = TextBox(ax_tb_xcenter, "", initial="0.0")

        ax_tb_ycenter = plt.axes([0.70, 0.13, 0.12, 0.03])
        tb_ycenter = TextBox(ax_tb_ycenter, "", initial="0.0")

        ax_tb_zoom = plt.axes([0.70, 0.07, 0.12, 0.03])
        tb_zoom = TextBox(ax_tb_zoom, "", initial="3.0")

        # Save button
        ax_save = plt.axes([0.85, 0.07, 0.10, 0.05])
        btn_save = Button(ax_save, "Save JPG")

        def update_plots() -> None:
            self.plot_2d(ax2d)
            self.plot_3d(ax3d)
            fig.canvas.draw_idle()

        def update_from_sliders(_) -> None:
            self.current_jacobi = slider_jacobi.val
            self.current_xcenter = slider_xcenter.val
            self.current_ycenter = slider_ycenter.val
            self.current_zoom = slider_zoom.val
            tb_jacobi.set_val(f"{self.current_jacobi:.10f}")
            tb_xcenter.set_val(f"{self.current_xcenter:.4f}")
            tb_ycenter.set_val(f"{self.current_ycenter:.4f}")
            tb_zoom.set_val(f"{self.current_zoom:.4f}")
            update_plots()

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
                update_plots()
            except ValueError:
                pass

        def update_ycenter_from_text(text: str) -> None:
            try:
                val = float(text)
                self.current_ycenter = val
                slider_ycenter.set_val(np.clip(val, -2.0, 2.0))
                update_plots()
            except ValueError:
                pass

        def update_zoom_from_text(text: str) -> None:
            try:
                val = float(text)
                if val > 0:
                    self.current_zoom = val
                    slider_zoom.set_val(np.clip(val, 0.001, 4.0))
                    update_plots()
            except ValueError:
                pass

        def save_image(_) -> None:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            default_name = f"zvc3d_C{self.current_jacobi:.10f}_x{self.current_xcenter:.2f}_y{self.current_ycenter:.2f}_z{self.current_zoom:.2f}.jpg"
            filepath = filedialog.asksaveasfilename(
                defaultextension=".jpg",
                filetypes=[("JPEG", "*.jpg"), ("PNG", "*.png"), ("All files", "*.*")],
                initialfile=default_name,
                title="Save ZVC Image"
            )
            root.destroy()
            if filepath:
                save_fig = plt.figure(figsize=(16, 9))
                save_ax2d = save_fig.add_subplot(121)
                save_ax3d = save_fig.add_subplot(122, projection="3d")
                self.plot_2d(save_ax2d)
                self.plot_3d(save_ax3d)
                save_fig.tight_layout()
                save_fig.savefig(filepath, dpi=300, bbox_inches="tight")
                plt.close(save_fig)
                print(f"Saved: {filepath}")

        slider_jacobi.on_changed(update_from_sliders)
        slider_xcenter.on_changed(update_from_sliders)
        slider_ycenter.on_changed(update_from_sliders)
        slider_zoom.on_changed(update_from_sliders)

        tb_jacobi.on_submit(update_jacobi_from_text)
        tb_xcenter.on_submit(update_xcenter_from_text)
        tb_ycenter.on_submit(update_ycenter_from_text)
        tb_zoom.on_submit(update_zoom_from_text)

        btn_save.on_clicked(save_image)

        # Initial plot
        update_plots()
        plt.show()


def main() -> None:
    """Main function."""
    plotter = ZVCPlotter3D(mu=MU_SUN_EARTH, resolution_2d=500, resolution_3d=80)
    plotter.run_interactive()


if __name__ == "__main__":
    main()

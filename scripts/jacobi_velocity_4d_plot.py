"""
4D Visualization of Asteroid Velocity Change and Jacobi Integral Relationship
in the Circular Restricted Three-Body Problem (CR3BP)

Axes configuration:
    - X-axis: Velocity change magnitude |Δv|
    - Y-axis: Jacobi integral change ΔC
    - Z-axis: Jacobi integral C
    - Color: Asteroid velocity |v| at that Jacobi integral

Features:
    - Interactive 3D rotation
    - Parameter sliders for mu, sample count, and perturbation range
    - Real-time data regeneration
    - View angle control
"""

from pathlib import Path
from typing import NamedTuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.widgets import Button, Slider
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


class AsteroidState(NamedTuple):
    """Data structure representing asteroid state (position and velocity)."""

    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float


def calc_r1(x: float, y: float, z: float, mu: float) -> float:
    """
    Calculate distance r1 from the third body to primary mass m1.

    Args:
        x: x-coordinate of the third body
        y: y-coordinate of the third body
        z: z-coordinate of the third body
        mu: Mass ratio mu = m2/(m1+m2)

    Returns:
        r1: Distance to primary mass m1
    """
    x1 = x + mu
    return np.sqrt(x1 * x1 + y * y + z * z)


def calc_r2(x: float, y: float, z: float, mu: float) -> float:
    """
    Calculate distance r2 from the third body to primary mass m2.

    Args:
        x: x-coordinate of the third body
        y: y-coordinate of the third body
        z: z-coordinate of the third body
        mu: Mass ratio mu = m2/(m1+m2)

    Returns:
        r2: Distance to primary mass m2
    """
    x2 = x - (1.0 - mu)
    return np.sqrt(x2 * x2 + y * y + z * z)


def calc_potential_u(x: float, y: float, z: float, mu: float) -> float:
    """
    Calculate effective potential U*.

    U* = 1/2(x^2 + y^2) + (1-mu)/r1 + mu/r2 + mu*(1-mu)/2

    Args:
        x: x-coordinate of the third body
        y: y-coordinate of the third body
        z: z-coordinate of the third body
        mu: Mass ratio mu = m2/(m1+m2)

    Returns:
        U*: Effective potential
    """
    r1 = calc_r1(x, y, z, mu)
    r2 = calc_r2(x, y, z, mu)
    if r1 == 0.0 or r2 == 0.0:
        raise ValueError("Position coincides with a primary.")
    return 0.5 * (x * x + y * y) + (1.0 - mu) / r1 + mu / r2 + mu * (1.0 - mu) * 0.5


def calc_jacobi_integral(state: AsteroidState, mu: float) -> float:
    """
    Calculate Jacobi integral C_J = 2*U* - v^2.

    Args:
        state: Asteroid state (position and velocity)
        mu: Mass ratio mu = m2/(m1+m2)

    Returns:
        C_J: Jacobi integral
    """
    v_sq = state.vx**2 + state.vy**2 + state.vz**2
    u_star = calc_potential_u(state.x, state.y, state.z, mu)
    return 2.0 * u_star - v_sq


def calc_velocity_magnitude(state: AsteroidState) -> float:
    """
    Calculate asteroid velocity magnitude.

    Args:
        state: Asteroid state

    Returns:
        |v|: Velocity magnitude
    """
    return np.sqrt(state.vx**2 + state.vy**2 + state.vz**2)


def generate_sample_data(
    n_samples: int = 1000,
    mu: float = 3.003e-6,
    perturb_range: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate sample data for visualization.

    Compute Jacobi integral changes when velocity perturbations are applied
    to asteroids at various positions and velocities.

    Args:
        n_samples: Number of samples to generate
        mu: Mass ratio (default: Sun-Earth system)
        perturb_range: Maximum velocity perturbation magnitude

    Returns:
        delta_v: Array of velocity change magnitudes |Δv|
        delta_c: Array of Jacobi integral changes ΔC
        jacobi_c: Array of Jacobi integrals C
        velocity: Array of asteroid velocities |v|
    """
    np.random.seed(42)

    delta_v_list = []
    delta_c_list = []
    jacobi_c_list = []
    velocity_list = []

    for _ in range(n_samples):
        # Generate random initial position (centered around L2)
        x = 1.01 + np.random.uniform(-0.05, 0.05)
        y = np.random.uniform(-0.1, 0.1)
        z = np.random.uniform(-0.02, 0.02)

        # Generate random initial velocity
        vx = np.random.uniform(-0.1, 0.1)
        vy = np.random.uniform(-0.1, 0.1)
        vz = np.random.uniform(-0.02, 0.02)

        initial_state = AsteroidState(x, y, z, vx, vy, vz)

        try:
            initial_jacobi = calc_jacobi_integral(initial_state, mu)
        except ValueError:
            continue

        # Apply random velocity perturbation
        dvx = np.random.uniform(-perturb_range, perturb_range)
        dvy = np.random.uniform(-perturb_range, perturb_range)
        dvz = np.random.uniform(-perturb_range * 0.2, perturb_range * 0.2)

        perturbed_state = AsteroidState(x, y, z, vx + dvx, vy + dvy, vz + dvz)

        try:
            perturbed_jacobi = calc_jacobi_integral(perturbed_state, mu)
            perturbed_velocity = calc_velocity_magnitude(perturbed_state)
        except ValueError:
            continue

        delta_v = np.sqrt(dvx**2 + dvy**2 + dvz**2)
        delta_c = perturbed_jacobi - initial_jacobi

        delta_v_list.append(delta_v)
        delta_c_list.append(delta_c)
        jacobi_c_list.append(perturbed_jacobi)
        velocity_list.append(perturbed_velocity)

    return (
        np.array(delta_v_list),
        np.array(delta_c_list),
        np.array(jacobi_c_list),
        np.array(velocity_list),
    )


class InteractivePlotter:
    """
    Interactive 4D plotter with GUI controls.

    Provides sliders for parameter adjustment and buttons for view control.
    """

    def __init__(self) -> None:
        """Initialize the interactive plotter."""
        self.fig = plt.figure(figsize=(14, 10))

        # Create 3D axis with adjusted position for controls
        self.ax = self.fig.add_axes([0.1, 0.25, 0.75, 0.7], projection="3d")

        # Initialize parameters
        self.mu = 3.003e-6
        self.n_samples = 1000
        self.perturb_range = 0.05
        self.elev = 20
        self.azim = 45

        # Generate initial data
        self.delta_v, self.delta_c, self.jacobi_c, self.velocity = generate_sample_data(
            self.n_samples, self.mu, self.perturb_range
        )

        # Store scatter plot reference
        self.scatter = None

        # Create GUI controls
        self._create_sliders()
        self._create_buttons()

        # Initial plot
        self._update_plot()

    def _create_sliders(self) -> None:
        """Create slider widgets for parameter control."""
        slider_color = "lightgoldenrodyellow"

        # Sample count slider
        ax_samples = self.fig.add_axes([0.15, 0.15, 0.55, 0.03], facecolor=slider_color)
        self.slider_samples = Slider(
            ax_samples,
            "Samples",
            100,
            5000,
            valinit=self.n_samples,
            valstep=100,
        )
        self.slider_samples.on_changed(self._on_samples_changed)

        # Perturbation range slider
        ax_perturb = self.fig.add_axes([0.15, 0.10, 0.55, 0.03], facecolor=slider_color)
        self.slider_perturb = Slider(
            ax_perturb,
            "Perturb Range",
            0.01,
            0.2,
            valinit=self.perturb_range,
        )
        self.slider_perturb.on_changed(self._on_perturb_changed)

        # Elevation angle slider
        ax_elev = self.fig.add_axes([0.15, 0.05, 0.25, 0.03], facecolor=slider_color)
        self.slider_elev = Slider(
            ax_elev,
            "Elevation",
            -90,
            90,
            valinit=self.elev,
        )
        self.slider_elev.on_changed(self._on_view_changed)

        # Azimuth angle slider
        ax_azim = self.fig.add_axes([0.45, 0.05, 0.25, 0.03], facecolor=slider_color)
        self.slider_azim = Slider(
            ax_azim,
            "Azimuth",
            -180,
            180,
            valinit=self.azim,
        )
        self.slider_azim.on_changed(self._on_view_changed)

    def _create_buttons(self) -> None:
        """Create button widgets for actions."""
        # Regenerate button
        ax_regen = self.fig.add_axes([0.75, 0.10, 0.12, 0.04])
        self.btn_regen = Button(ax_regen, "Regenerate", color="lightblue", hovercolor="skyblue")
        self.btn_regen.on_clicked(self._on_regenerate)

        # Reset view button
        ax_reset = self.fig.add_axes([0.75, 0.05, 0.12, 0.04])
        self.btn_reset = Button(ax_reset, "Reset View", color="lightgreen", hovercolor="palegreen")
        self.btn_reset.on_clicked(self._on_reset_view)

        # Save button
        ax_save = self.fig.add_axes([0.88, 0.10, 0.1, 0.04])
        self.btn_save = Button(ax_save, "Save", color="lightyellow", hovercolor="khaki")
        self.btn_save.on_clicked(self._on_save)

    def _update_plot(self) -> None:
        """Update the 3D scatter plot with current data."""
        self.ax.clear()

        # Set up colormap normalization
        if len(self.velocity) > 0:
            norm = Normalize(vmin=self.velocity.min(), vmax=self.velocity.max())
        else:
            norm = Normalize(vmin=0, vmax=1)

        # Create scatter plot
        self.scatter = self.ax.scatter(
            self.delta_v,
            self.delta_c,
            self.jacobi_c,
            c=self.velocity,
            cmap=cm.viridis,
            norm=norm,
            s=15,
            alpha=0.7,
            edgecolors="none",
        )

        # Set axis labels
        self.ax.set_xlabel(r"Velocity Change $|\Delta v|$ [-]", fontsize=11, labelpad=8)
        self.ax.set_ylabel(r"Jacobi Integral Change $\Delta C$ [-]", fontsize=11, labelpad=8)
        self.ax.set_zlabel(r"Jacobi Integral $C$ [-]", fontsize=11, labelpad=8)

        # Set title
        self.ax.set_title(
            f"CR3BP: Velocity Change vs Jacobi Integral (n={len(self.delta_v)})\n"
            f"Color: Asteroid Velocity $|v|$",
            fontsize=12,
        )

        # Set view angle
        self.ax.view_init(elev=self.elev, azim=self.azim)

        # Enable grid
        self.ax.grid(True, alpha=0.3)

        # Update colorbar (remove old one if exists)
        if hasattr(self, "cbar") and self.cbar is not None:
            self.cbar.remove()

        if len(self.velocity) > 0:
            self.cbar = self.fig.colorbar(
                self.scatter, ax=self.ax, shrink=0.5, aspect=15, pad=0.12
            )
            self.cbar.set_label(r"Asteroid Velocity $|v|$ [-]", fontsize=10)
        else:
            self.cbar = None

        self.fig.canvas.draw_idle()

    def _on_samples_changed(self, val: float) -> None:
        """Handle sample count slider change."""
        self.n_samples = int(val)

    def _on_perturb_changed(self, val: float) -> None:
        """Handle perturbation range slider change."""
        self.perturb_range = val

    def _on_view_changed(self, val: float) -> None:
        """Handle view angle slider change."""
        self.elev = self.slider_elev.val
        self.azim = self.slider_azim.val
        self.ax.view_init(elev=self.elev, azim=self.azim)
        self.fig.canvas.draw_idle()

    def _on_regenerate(self, event) -> None:
        """Handle regenerate button click."""
        print(f"Regenerating with n={self.n_samples}, perturb={self.perturb_range:.3f}...")
        self.delta_v, self.delta_c, self.jacobi_c, self.velocity = generate_sample_data(
            self.n_samples, self.mu, self.perturb_range
        )
        self._update_plot()
        print(f"Generated {len(self.delta_v)} data points.")

    def _on_reset_view(self, event) -> None:
        """Handle reset view button click."""
        self.slider_elev.set_val(20)
        self.slider_azim.set_val(45)

    def _on_save(self, event) -> None:
        """Handle save button click."""
        output_dir = Path(__file__).parent / "output"
        output_dir.mkdir(exist_ok=True)
        output_path = output_dir / "jacobi_velocity_4d_plot.png"
        self.fig.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Saved to: {output_path}")

    def show(self) -> None:
        """Display the interactive plot."""
        plt.show()


def main() -> None:
    """Main function to launch the interactive plotter."""
    print("Launching interactive 4D plotter...")
    print("Controls:")
    print("  - Samples: Adjust number of data points")
    print("  - Perturb Range: Adjust velocity perturbation magnitude")
    print("  - Elevation/Azimuth: Adjust view angle")
    print("  - Regenerate: Generate new random data")
    print("  - Reset View: Reset to default view angle")
    print("  - Save: Save current plot to file")
    print("  - Mouse drag: Rotate 3D view")

    plotter = InteractivePlotter()
    plotter.show()


if __name__ == "__main__":
    main()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def potential(x, y, mu):
    r1 = np.sqrt((x + mu) ** 2 + y**2)
    r2 = np.sqrt((x - 1 + mu) ** 2 + y**2)
    return x**2 + y**2 + 2 * (1 - mu) / r1 + 2 * mu / r2


def zero_velocity_curves(mu, C_values):
    x = np.linspace(0.97, 1.03, 500)
    y = np.linspace(-0.03, 0.03, 500)
    X, Y = np.meshgrid(x, y)

    Z = potential(X, Y, mu)

    plt.figure(figsize=(12, 10))
    colors = ["#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"]
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list("custom", colors, N=n_bins)

    contour = plt.contourf(X, Y, Z, levels=C_values, cmap=cmap)
    plt.colorbar(contour, label="Jacobi Integral")

    plt.plot([0.97, 1.03], [0, 0], "k--")
    plt.title(f"Zero Velocity Curves (μ = {mu})")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.grid(True)
    plt.show()


# パラメータの設定
mu = 0.000003003  # 質量比
C_values = np.linspace(3.00, 3.0011, 100)  # ヤコビ積分の値の配列
zero_velocity_curves(mu, C_values)

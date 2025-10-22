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

    plt.figure(figsize=(14, 14))

    plt.contourf(X, Y, Z, levels=[-np.inf, C_values], colors="blue", alpha=0.5)

    plt.plot([0.97, 1.03], [0, 0], "k--")
    plt.title(f"Zero Velocity Curves (C_j = {C_values})")
    plt.xlabel("x(au)")
    plt.ylabel("y(au)")
    plt.axis("equal")
    plt.grid(True)
    plt.show()


# パラメータの設定
mu = 0.000003003  # 質量比
# C_values = np.arange(3.0009, 3.001, 0.00001)  # ヤコビ定数の値
C_values = 3.00089055

# C_values = 3.000891
zero_velocity_curves(mu, C_values)

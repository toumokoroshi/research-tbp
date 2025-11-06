import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# def potential(x, y, mu):
#     r = np.sqrt(x**2 + y**2)
#     return -1 / r


def potential(x, y, mu):
    r1 = np.sqrt((x + mu) ** 2 + y**2)
    r2 = np.sqrt((x - 1 + mu) ** 2 + y**2)
    return x**2 + y**2 + 2 * (1 - mu) / r1 + 2 * mu / r2


def zero_velocity_curves(mu, C_values):
    # x = np.linspace(0.97, 1.03, 500)
    # y = np.linspace(-0.03, 0.03, 500)
    # x = np.linspace(-1.0, 1.03, 500)
    x = np.linspace(-1.5, 1.5, 1000)
    y = np.linspace(-1.5, 1.5, 1000)
    X, Y = np.meshgrid(x, y)

    Z = potential(X, Y, mu)

    plt.figure(figsize=(14, 14))

    # plt.contourf(X, Y, Z, levels=[C_values, np.inf], colors="blue", alpha=0.5)
    plt.contourf(X, Y, Z, levels=[-np.inf, C_values], colors="blue", alpha=0.5)

    # plt.plot([0.97, 1.03], [0, 0], "k--")
    # plt.plot([-1.0, 1.03], [0, 0], "k--")
    plt.title(f"Zero Velocity Curves (C_j = {C_values})")
    plt.xlabel("x(au)")
    plt.ylabel("y(au)")
    plt.axis("equal")
    plt.grid(True)
    plt.show()


# パラメータの設定
gm_sun = 1.32712440041279419e20  # // heliocentric gravitational parameter (m3 s-2)
gm_earth = 3.98600435507e14  # // geocentric gravitational parameter (m3 s-2)
mu_sun_earth = gm_earth / (gm_sun + gm_earth)
mu = 0  # 質量比
# mu = 0.000003003  # 質量比
# C_values = np.arange(3.0009, 3.001, 0.00001)  # ヤコビ定数の値
# C_values = 3.0

C_values = 4
zero_velocity_curves(mu_sun_earth, C_values)

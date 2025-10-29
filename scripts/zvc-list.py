import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def effective_potential(x, y, mu):
    """
    制限三体問題の有効ポテンシャル
    mu: 質量比パラメータ（小さい天体の質量 / 全質量）
    """
    # 二つの主天体の位置
    r1 = np.sqrt((x + mu) ** 2 + y**2)
    r2 = np.sqrt((x - 1 + mu) ** 2 + y**2)

    # ゼロ除算を避ける
    r1 = np.maximum(r1, 1e-10)
    r2 = np.maximum(r2, 1e-10)

    # 有効ポテンシャル
    U = 0.5 * (x**2 + y**2) + (1 - mu) / r1 + mu / r2

    return U


def plot_zero_velocity_curves(mu=0.01, jacobi_constants=None):
    """
    ゼロ速度曲線を描画（天体2周辺のズーム付き）
    mu: 質量比（例: 地球-月系なら約0.012）
    jacobi_constants: ヤコビ定数のリスト
    """
    # グリッドの作成（全体図用）
    x = np.linspace(-1.5, 1.5, 500)
    y = np.linspace(-1.5, 1.5, 500)
    X, Y = np.meshgrid(x, y)

    # 有効ポテンシャルの計算
    U = effective_potential(X, Y, mu)

    # ズーム領域用のグリッド（天体2周辺）
    zoom_range = 0.05  # ズーム範囲
    x_zoom = np.linspace(1 - mu - zoom_range, 1 - mu + zoom_range, 300)
    y_zoom = np.linspace(-zoom_range, zoom_range, 300)
    X_zoom, Y_zoom = np.meshgrid(x_zoom, y_zoom)
    U_zoom = effective_potential(X_zoom, Y_zoom, mu)

    # ヤコビ定数の設定（デフォルト値）
    if jacobi_constants is None:
        jacobi_constants = [3.0, 3.2, 3.4, 3.6, 3.8]

    # プロットの作成
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for idx, C in enumerate(jacobi_constants):
        if idx >= 6:
            break

        ax = axes[idx]

        # ゼロ速度曲線: 2U = C の等高線
        # 運動可能領域: 2U <= C (速度^2 = C - 2U >= 0)
        # 運動不可能領域: 2U > C

        forbidden_region = 2 * U > C
        forbidden_region_zoom = 2 * U_zoom > C

        # 運動不可能な領域を塗りつぶし（全体図）
        ax.contourf(
            X,
            Y,
            forbidden_region.astype(int),
            levels=[-0.5, 0.5],
            colors=["lightcoral"],
            alpha=0.6,
        )

        # ゼロ速度曲線を描画（全体図）
        contours = ax.contour(X, Y, 2 * U, levels=[C], colors="red", linewidths=2)

        # 二つの主天体の位置をプロット
        ax.plot(-mu, 0, "o", color="orange", markersize=10, label="Sun")
        ax.plot(1 - mu, 0, "o", color="skyblue", markersize=3, label="Earth")

        # # ラグランジュ点L1, L2, L3の概略位置（簡易計算）
        # L1_x = 1 - mu - (mu / 3) ** (1 / 3)
        # L2_x = 1 - mu + (mu / 3) ** (1 / 3)
        # L3_x = -mu - (1 - mu) * (7 * mu / 12)

        # ax.plot(
        #     [L1_x, L2_x, L3_x],
        #     [0, 0, 0],
        #     "x",
        #     color="orange",
        #     markersize=8,
        #     label="L1,L2,L3",
        # )

        # ズーム領域を示す矩形を描画
        from matplotlib.patches import Rectangle, FancyBboxPatch, ConnectionPatch

        zoom_rect = Rectangle(
            (1 - mu - zoom_range, -zoom_range),
            2 * zoom_range,
            2 * zoom_range,
            linewidth=2,
            edgecolor="darkblue",
            facecolor="none",
            linestyle="-",
        )
        ax.add_patch(zoom_rect)

        # インセット軸（ズーム図）を作成 - 軸の外側に配置
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        ax_zoom = inset_axes(
            ax,
            width="50%",
            height="50%",
            bbox_to_anchor=(0.7, 0.5, 0.9, 0.9),
            bbox_transform=ax.transAxes,
            loc="center left",
            borderpad=0,
        )

        # ズーム図に背景色を設定
        ax_zoom.set_facecolor("#FFFFF0")

        # ズーム図に運動不可能領域を塗りつぶし
        ax_zoom.contourf(
            X_zoom,
            Y_zoom,
            forbidden_region_zoom.astype(int),
            levels=[-0.5, 0.5],
            colors=["lightcoral"],
            alpha=0.6,
        )

        # ズーム図にゼロ速度曲線を描画
        ax_zoom.contour(
            X_zoom, Y_zoom, 2 * U_zoom, levels=[C], colors="red", linewidths=2
        )

        # ズーム図に天体2とラグランジュ点を表示
        ax_zoom.plot(1 - mu, 0, "o", color="skyblue", markersize=3)

        ax_zoom.set_xlim(1 - mu - zoom_range, 1 - mu + zoom_range)
        ax_zoom.set_ylim(-zoom_range, zoom_range)
        ax_zoom.set_aspect("equal")
        ax_zoom.grid(True, alpha=0.3, linewidth=0.5)
        ax_zoom.tick_params(labelsize=8)

        # ズーム図の枠を強調
        for spine in ax_zoom.spines.values():
            spine.set_edgecolor("darkblue")
            spine.set_linewidth(2.5)

        # 元の領域とズーム図を線で繋ぐ
        # 右下の角を左下に繋ぐ
        con1 = ConnectionPatch(
            (1 - mu - zoom_range, -zoom_range),
            (0, 0),
            "data",
            "axes fraction",
            axesA=ax,
            axesB=ax_zoom,
            color="darkblue",
            linewidth=1.5,
            linestyle="-",
            alpha=0.7,
        )
        ax.add_artist(con1)

        # 右上の角を左上に繋ぐ
        con2 = ConnectionPatch(
            (1 - mu + zoom_range, -zoom_range),
            (1, 0),
            "data",
            "axes fraction",
            axesA=ax,
            axesB=ax_zoom,
            color="darkblue",
            linewidth=1.5,
            linestyle="-",
            alpha=0.7,
        )
        ax.add_artist(con2)

        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(f"C_j = {C:.6f}")
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    # 最後のサブプロットを非表示
    if len(jacobi_constants) < 6:
        axes[5].axis("off")

    plt.suptitle(
        f"ZVC and forbidden regions \n" + f"μ = {mu:e}",
        fontsize=14,
        y=0.98,
    )

    plt.tight_layout()
    plt.show()


# 実行例
print()

# 地球-月系の質量比
mu_earth_moon = 0.01215
gm_sun = 1.32712440041279419e20  # // heliocentric gravitational parameter (m3 s-2)
gm_earth = 3.98600435507e14  # // geocentric gravitational parameter (m3 s-2)
mu_sun_earth = gm_earth / (gm_sun + gm_earth)
# 様々なヤコビ定数でプロット
jacobi_values = [3.0, 3.0002, 3.0004, 3.0006, 3.0008, 3.001]
plot_zero_velocity_curves(mu=mu_sun_earth, jacobi_constants=jacobi_values)

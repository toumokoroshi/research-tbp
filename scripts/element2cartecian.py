import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plot_hill_sphere_au():
    # --- 定数設定 (無次元化) ---
    # 距離の基準: 1.0 = 1 AU (太陽-地球間距離)
    a_AU = 1.0

    # 質量の基準
    M_sun = 1.989e30    # 太陽質量 (kg)
    M_earth = 5.972e24  # 地球質量 (kg)

    # ヒル球の半径計算 (単位: AU)
    # 公式: r = a * (m / 3M)^(1/3)
    r_hill_au = a_AU * (M_earth / (3 * M_sun))**(1/3)

    print(f"ヒル球の半径: {r_hill_au:.5f} AU")

    # --- 参考データ: 月の軌道 ---
    # 月までの距離: 約384,400 km
    # 1 AU: 約149,600,000 km
    km_per_au = 1.496e8
    r_moon_au = 384400 / km_per_au

    # --- プロット作成 ---
    fig, ax = plt.subplots(figsize=(8, 8))

    # 1. ヒル球 (Hill Sphere)
    hill_circle = patches.Circle((0, 0), r_hill_au,
                                 color='green', alpha=0.1,
                                 label='ヒル球 (Hill Sphere)')
    ax.add_patch(hill_circle)
    ax.add_patch(patches.Circle((0, 0), r_hill_au, color='green', fill=False, linestyle='--'))

    # 2. 月の軌道 (Moon Orbit)
    moon_circle = patches.Circle((0, 0), r_moon_au,
                                 edgecolor='gray', linestyle=':', fill=False,
                                 label='月の軌道 (Moon Orbit)')
    ax.add_patch(moon_circle)

    # 3. 地球 (Earth)
    ax.plot(0, 0, marker='o', color='blue', markersize=8, label='地球 (Earth)')

    # 4. ラグランジュ点 L1, L2 (近似)
    ax.scatter([-r_hill_au, r_hill_au], [0, 0], color='red', marker='x', label='L1 / L2 (近似)')

    # 5. 太陽方向の矢印
    ax.annotate('太陽の方向 (-1.0 AU)', xy=(-r_hill_au, 0), xytext=(-r_hill_au*1.4, 0),
                arrowprops=dict(facecolor='orange', arrowstyle='->', color='orange'),
                ha='right', va='center', color='darkorange', weight='bold')

    # --- グラフの装飾 ---
    limit = r_hill_au * 1.5
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect('equal')

    ax.set_xlabel('地球からの距離 (AU)')
    ax.set_ylabel('地球からの距離 (AU)')
    ax.set_title(f'地球のヒル球 (無次元化距離)\n半径 $r_{{Hill}} \\approx {r_hill_au:.5f}$ AU', fontsize=14)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_hill_sphere_au()

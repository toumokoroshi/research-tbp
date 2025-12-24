import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import ConnectionPatch
from scipy.integrate import solve_ivp

# --- パラメータ設定 (太陽-地球系) ---
mu = 3.003e-6  # 太陽-地球の質量比
C_target = 3.0006  # ヤコビ定数


# --- 運動方程式 (CR3BP) ---
def equations(t, state, mu):
    x, y, vx, vy = state
    r1 = ((x + mu) ** 2 + y**2) ** 1.5
    r2 = ((x - (1 - mu)) ** 2 + y**2) ** 1.5

    ax = 2 * vy + x - (1 - mu) * (x + mu) / r1 - mu * (x - (1 - mu)) / r2
    ay = -2 * vx + y - (1 - mu) * y / r1 - mu * y / r2
    return [vx, vy, ax, ay]


# --- イベント: y=0 面通過 ---
def event_poincare(t, state, mu):
    return state[1]


event_poincare.direction = 1  # 負から正へ (y=0を上に抜ける時)
event_poincare.terminal = False


# --- イベント: 領域外への脱出（計算打ち切り用） ---
def event_escape(t, state, mu):
    x, y, _, _ = state
    r2 = np.sqrt((x - (1 - mu)) ** 2 + y**2)
    if r2 < 1e-12:
        return -1  # 地球衝突
    return 1.0


event_escape.terminal = True


# --- 初期速度 vy の計算 ---
def get_vy(x, C, mu):
    y = 0
    r1 = np.sqrt((x + mu) ** 2 + y**2)
    r2 = np.sqrt((x - (1 - mu)) ** 2 + y**2)
    Omega = 0.5 * (x**2 + y**2) + (1 - mu) / r1 + mu / r2
    v2 = 2 * Omega - C
    if v2 < 0:
        return None
    return np.sqrt(v2)


# --- シミュレーション実行関数 ---
def run_simulation(x0, t_max=300):
    """指定されたx0でシミュレーションを実行し、解とイベントを返す"""
    vy0 = get_vy(x0, C_target, mu)
    if vy0 is None:
        return None, None

    sol = solve_ivp(
        equations,
        [0, t_max],
        [x0, 0, 0, vy0],
        args=(mu,),
        events=[event_poincare, event_escape],
        rtol=1e-11,
        atol=1e-11,
        max_step=0.01,
    )
    return sol, vy0


# 計算範囲の設定
x_values = np.linspace(0.97, 1.03, 120)  # 点数を少し増やしても良いでしょう

# 結果格納用のリスト
poincare_results = []

print(f"Generating Poincare Map for C={C_target}...")
print(f"Total points to scan: {len(x_values)}")

for i, x0 in enumerate(x_values):
    sol, _ = run_simulation(x0, t_max=300)

    if sol is not None and len(sol.t_events[0]) > 0:
        # イベント（断面通過点）の座標だけ保存する
        # ys: [x, y, vx, vy] の配列
        ys = sol.y_events[0]
        poincare_results.append(ys)

    if i % 20 == 0:
        print(f"Progress: {i}/{len(x_values)}")

print("Calculation Done.")

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import ConnectionPatch
from matplotlib.ticker import MaxNLocator  # <--- ★目盛りを間引くための機能を追加

# --- 注目したい軌道の再設定 ---
x_target_periodic = 0.996
x_target_chaotic = 0.9992278055234482

# --- フォント設定 ---
plt.rcParams.update(
    {
        "font.size": 18,
        "axes.titlesize": 16,
        "axes.labelsize": 16,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "figure.titlesize": 20,
    }
)

# --- キャンバス準備 ---
fig = plt.figure(figsize=(14, 8))
gs = GridSpec(2, 2, width_ratios=[1.5, 1])

ax_map = fig.add_subplot(gs[:, 0])  # 左: ポアンカレマップ
ax_orbit1 = fig.add_subplot(gs[0, 1])  # 右上: ターゲット軌道1
ax_orbit2 = fig.add_subplot(gs[1, 1])  # 右下: ターゲット軌道2

# 1. ポアンカレマップの描画
print("Plotting Map...")
for ys in poincare_results:
    ax_map.plot(ys[:, 0], ys[:, 2], ".", markersize=1, color="tab:blue", alpha=0.4)

# 2. 特定軌道の個別計算と描画
targets = [
    (
        x_target_periodic,
        "Regular",
        ax_orbit1,
        f"Regular Orbit (x0={x_target_periodic})",
    ),
    (x_target_chaotic, "Chaotic", ax_orbit2, f"Chaotic Orbit (x0={x_target_chaotic})"),
]

for x_tgt, label, ax_orb, title in targets:
    sol, vy0 = run_simulation(x_tgt, t_max=200)

    if sol is not None:
        # A. マップ上に強調点
        ax_map.plot(
            x_tgt, 0, "o", color="red", markersize=6, zorder=10, markeredgecolor="black"
        )

        # B. 実軌道プロット
        ax_orb.plot(sol.y[0], sol.y[1], linewidth=1.0, color="tab:orange")

        # --- ★文字位置調整エリア (右側のグラフ) ---
        ax_orb.set_title(title, pad=10)  # pad: タイトルとグラフの距離
        ax_orb.set_xlabel("x", labelpad=5)  # labelpad: 軸ラベルと目盛りの距離
        ax_orb.set_ylabel("y", labelpad=5, rotation=0)  # rotation: 文字の回転

        ax_orb.axis("equal")
        ax_orb.grid(True, alpha=0.3)

        # --- ★目盛りの間引き (右側のグラフ) ---
        # nbins=4 は「目盛りを最大4つまでにする」という意味です
        ax_orb.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax_orb.yaxis.set_major_locator(MaxNLocator(nbins=4))

        # C. 矢印による接続
        con = ConnectionPatch(
            xyA=(x_tgt, 0),
            coordsA=ax_map.transData,
            xyB=(0, 0.5),
            coordsB=ax_orb.transAxes,
            arrowstyle="-|>",
            color="black",
            alpha=0.8,
            lw=1.5,
        )
        fig.add_artist(con)

# --- ポアンカレマップの装飾 ---
# --- ★文字位置調整エリア (左側のグラフ) ---
ax_map.set_title(
    f"Poincaré Section (Sun-Earth L1), C={C_target}", y=1.02
)  # y: タイトルの縦位置(1.0が枠の上端)
ax_map.set_xlabel("x (ND)", labelpad=10)  # labelpad: 数値を大きくすると下に下がる
ax_map.set_ylabel("vx (ND)", labelpad=10)

ax_map.set_xlim(0.985, 1.005)
ax_map.set_ylim(-0.02, 0.02)
ax_map.grid(True, alpha=0.3)

# --- ★目盛りの間引き (左側のグラフ) ---
ax_map.xaxis.set_major_locator(MaxNLocator(nbins=5))  # x軸の目盛りを5つ程度に制限
ax_map.yaxis.set_major_locator(MaxNLocator(nbins=6))  # y軸の目盛りを6つ程度に制限

# --- ★グラフ同士の距離を近づける設定 ---
# tight_layoutをした後に、subplots_adjustで微調整します
plt.tight_layout()
plt.subplots_adjust(wspace=0.1, hspace=0.2)
# wspace: 横の隙間 (デフォルトは0.2くらい)
# hspace: 縦の隙間 (デフォルトは0.2くらい)

plt.show()

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint


# --- 1. 物理モデルの定義 (無次元化) ---
def double_pendulum_derivs_dimensionless(state, tau, l, mu):
    """
    無次元化された二重振り子の運動方程式 (Hamiltonian系に由来)

    Args:
        state: [theta1, omega1, theta2, omega2]
        tau: 無次元時間
        l: 長さの比 (L2 / L1)
        mu: 質量の比 (m2 / m1)

    Returns:
        状態ベクトルの時間微分
    """
    th1, w1, th2, w2 = state
    delta = th2 - th1
    cos_delta = np.cos(delta)
    sin_delta = np.sin(delta)

    den1 = (1 + mu) - mu * cos_delta**2
    den2 = l * den1

    dydt = np.zeros_like(state)
    dydt[0] = w1
    dydt[2] = w2

    num1 = (
        mu * w1**2 * sin_delta * cos_delta
        + mu * np.sin(th2) * cos_delta
        + mu * l * w2**2 * sin_delta
        - (1 + mu) * np.sin(th1)
    )
    dydt[1] = num1 / den1

    num2 = -mu * l * w2**2 * sin_delta * cos_delta + (1 + mu) * (
        np.sin(th1) * cos_delta - w1**2 * sin_delta - np.sin(th2)
    )
    dydt[3] = num2 / den2

    return dydt


# --- 2. 無次元パラメータ設定 ---
l = 1.0  # 長さの比 (L2 / L1)
mu = 1.0  # 質量の比 (m2 / m1)

# 時間設定
dt = 0.04  # 無次元時間刻み -> アニメーションのフレームレートに影響
t_max = 20.0  # シミュレーション終了時間
t = np.arange(0, t_max, dt)

# --- 3. 初期条件の設定 (2つのケース) ---
# ケース1: 基準となる初期状態
# [theta1, omega1, theta2, omega2] (角度はラジアン)
th1 = 180
th2 = 90
state1_init = np.radians([th1, 0, th2, 0])

# ケース2: ケース1から「わずかに」ずらした状態
# theta1を 0.001度 だけ変える
state2_init = np.radians([th1 + 0.001, 0, th2, 0])

# --- 4. 運動方程式を解く ---
print("シミュレーション計算中...")
y1 = odeint(double_pendulum_derivs_dimensionless, state1_init, t, args=(l, mu))
y2 = odeint(double_pendulum_derivs_dimensionless, state2_init, t, args=(l, mu))


# --- 5. 座標変換 (極座標 -> 直交座標) ---
# 無次元化: L1 = 1 (基準長), L2 = l * L1 = l
def get_coords(y, l):
    theta1 = y[:, 0]
    theta2 = y[:, 2]

    x1 = np.sin(theta1)
    y1 = -np.cos(theta1)

    x2 = x1 + l * np.sin(theta2)
    y2 = y1 - l * np.cos(theta2)
    return x1, y1, x2, y2


x1_a, y1_a, x2_a, y2_a = get_coords(y1, l)
x1_b, y1_b, x2_b, y2_b = get_coords(y2, l)

# --- 6. アニメーション設定 ---
fig, ax = plt.subplots(figsize=(8, 8))
lim = 1 + l + 0.2
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect("equal")
ax.grid()
ax.set_title(f"Double Pendulum (Dimensionless: l={l}, μ={mu})")

# 描画オブジェクトの初期化
# 振り子1 (Blue)
(line_a,) = ax.plot([], [], "o-", lw=2, color="blue", label="Initial A")
(trace_a,) = ax.plot([], [], "-", lw=1, alpha=0.5, color="blue")  # 軌跡

# 振り子2 (Red) - 重ねて表示
(line_b,) = ax.plot([], [], "o-", lw=2, color="red", label="Initial B (tiny offset)")
(trace_b,) = ax.plot([], [], "-", lw=1, alpha=0.5, color="red")  # 軌跡

# 時間表示用テキスト
time_template = "Time = {:.1f}s"
time_text = ax.text(0.05, 0.9, "", transform=ax.transAxes)

ax.legend(loc="upper right")


def init():
    line_a.set_data([], [])
    trace_a.set_data([], [])
    line_b.set_data([], [])
    trace_b.set_data([], [])
    time_text.set_text("")
    return line_a, trace_a, line_b, trace_b, time_text


def update(i):
    # --- 振り子Aの更新 ---
    thisx_a = [0, x1_a[i], x2_a[i]]
    thisy_a = [0, y1_a[i], y2_a[i]]
    line_a.set_data(thisx_a, thisy_a)

    # 軌跡A (過去の履歴すべてを表示、長すぎる場合は [-100:] などでスライス)
    trace_a.set_data(x2_a[:i], y2_a[:i])

    # --- 振り子Bの更新 ---
    thisx_b = [0, x1_b[i], x2_b[i]]
    thisy_b = [0, y1_b[i], y2_b[i]]
    line_b.set_data(thisx_b, thisy_b)

    # 軌跡B
    trace_b.set_data(x2_b[:i], y2_b[:i])

    time_text.set_text(time_template.format(i * dt))
    return line_a, trace_a, line_b, trace_b, time_text


# アニメーション作成
ani = animation.FuncAnimation(
    fig, update, frames=len(t), interval=dt * 1000, blit=True, init_func=init
)

plt.show()

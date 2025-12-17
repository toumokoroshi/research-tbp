import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint


# --- 1. 物理モデルの定義 ---
def derivs(state, t, L1, L2, m1, m2, g):
    """
    二重振り子の運動方程式 (ラグランジュ力学より導出)
    state: [theta1, omega1, theta2, omega2]
    """
    theta1, w1, theta2, w2 = state

    delta = theta2 - theta1
    den1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta) * np.cos(delta)
    den2 = (L2 / L1) * den1

    # dw1/dt (角加速度1)
    dw1 = (
        m2 * L1 * w1 * w1 * np.sin(delta) * np.cos(delta)
        + m2 * g * np.sin(theta2) * np.cos(delta)
        + m2 * L2 * w2 * w2 * np.sin(delta)
        - (m1 + m2) * g * np.sin(theta1)
    ) / den1

    # dw2/dt (角加速度2)
    dw2 = (
        -m2 * L2 * w2 * w2 * np.sin(delta) * np.cos(delta)
        + (m1 + m2) * g * np.sin(theta1) * np.cos(delta)
        - (m1 + m2) * L1 * w1 * w1 * np.sin(delta)
        - (m1 + m2) * g * np.sin(theta2)
    ) / den2

    return [w1, dw1, w2, dw2]


# --- 2. パラメータ設定 ---
# 物理定数
L1, L2 = 1.0, 1.0  # ロッドの長さ (m)
m1, m2 = 1.0, 1.0  # 重りの質量 (kg)
g = 9.8  # 重力加速度 (m/s^2)

# 時間設定
dt = 0.04  # 時間刻み (秒) -> アニメーションのフレームレートに影響
t_max = 20.0  # シミュレーション終了時間
t = np.arange(0, t_max, dt)

# --- 3. 初期条件の設定 (2つのケース) ---
# ケース1: 基準となる初期状態
# [theta1, omega1, theta2, omega2] (角度はラジアン, 120度からスタート)
th1 = 180
th2 = 90
state1_init = np.radians([th1, 0, th2, 0])

# ケース2: ケース1から「わずかに」ずらした状態
# theta1を 0.1度 だけ変える
state2_init = np.radians([th1 + 0.001, 0, th2, 0])

# --- 4. 運動方程式を解く ---
print("シミュレーション計算中...")
# args=(L1, L2, m1, m2, g) で定数を渡す
y1 = odeint(derivs, state1_init, t, args=(L1, L2, m1, m2, g))
y2 = odeint(derivs, state2_init, t, args=(L1, L2, m1, m2, g))


# --- 5. 座標変換 (極座標 -> 直交座標) ---
def get_coords(y, L1, L2):
    theta1 = y[:, 0]
    theta2 = y[:, 2]

    x1 = L1 * np.sin(theta1)
    y1 = -L1 * np.cos(theta1)

    x2 = x1 + L2 * np.sin(theta2)
    y2 = y1 - L2 * np.cos(theta2)
    return x1, y1, x2, y2


x1_a, y1_a, x2_a, y2_a = get_coords(y1, L1, L2)
x1_b, y1_b, x2_b, y2_b = get_coords(y2, L1, L2)

# --- 6. アニメーション設定 ---
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-(L1 + L2 + 0.2), (L1 + L2 + 0.2))
ax.set_ylim(-(L1 + L2 + 0.2), (L1 + L2 + 0.2))
ax.set_aspect("equal")
ax.grid()
ax.set_title("Double Pendulum Simulation (Chaos)")

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

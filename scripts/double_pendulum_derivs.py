import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint

# --- パラメータ設定 ---
g = 9.81
L1, L2 = 1.0, 1.0  # ロッドの長さ [m]
m1, m2 = 1.0, 1.0  # おもりの質量 [kg]
t_max = 20.0  # シミュレーション時間 [s]
dt = 0.05  # 時間刻み

# 初期状態: [theta1, omega1, theta2, omega2]
init_state = [np.pi / 2 + 0.01, 0.0, np.pi, 0.0]


# --- 運動方程式 (二重振り子) ---
def double_pendulum_derivs(state, t, L1, L2, m1, m2, g):
    th1, w1, th2, w2 = state

    delta = th2 - th1
    den1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta) * np.cos(delta)
    den2 = (L2 / L1) * den1

    dydt = np.zeros_like(state)
    dydt[0] = w1
    dydt[2] = w2

    num1 = (
        m2 * L1 * w1**2 * np.sin(delta) * np.cos(delta)
        + m2 * g * np.sin(th2) * np.cos(delta)
        + m2 * L2 * w2**2 * np.sin(delta)
        - (m1 + m2) * g * np.sin(th1)
    )
    dydt[1] = num1 / den1

    num2 = -m2 * L2 * w2**2 * np.sin(delta) * np.cos(delta) + (m1 + m2) * (
        g * np.sin(th1) * np.cos(delta) - L1 * w1**2 * np.sin(delta) - g * np.sin(th2)
    )
    dydt[3] = num2 / den2

    return dydt


t = np.arange(0, t_max, dt)
y = odeint(double_pendulum_derivs, init_state, t, args=(L1, L2, m1, m2, g))

th1 = y[:, 0]
w1 = y[:, 1]
th2 = y[:, 2]
w2 = y[:, 3]

# 直交座標へ変換
x1 = L1 * np.sin(th1)
y1 = -L1 * np.cos(th1)
x2 = x1 + L2 * np.sin(th2)
y2 = y1 - L2 * np.cos(th2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

ax1.set_title("Real Space Motion")
ax1.set_xlabel("x [m]")
ax1.set_ylabel("y [m]")
ax1.set_xlim(-(L1 + L2) * 1.1, (L1 + L2) * 1.1)
ax1.set_ylim(-(L1 + L2) * 1.1, (L1 + L2) * 1.1)
ax1.set_aspect("equal")
ax1.grid()

(line,) = ax1.plot([], [], "o-", lw=2, color="black")  # ロッド
(mass1,) = ax1.plot([], [], "o", color="blue", markersize=8, label="Mass 1")
(mass2,) = ax1.plot([], [], "o", color="red", markersize=8, label="Mass 2")
(trace,) = ax1.plot([], [], "-", lw=0.5, color="gray", alpha=0.5)  # 軌跡
time_text = ax1.text(0.05, 0.9, "", transform=ax1.transAxes)
ax1.legend(loc="lower left")

ax2.set_title("Phase Space Projections")
ax2.set_xlabel(r"Angle $\theta$ [rad]")
ax2.set_ylabel(r"Angular Velocity $\omega$ [rad/s]")
ax2.set_xlim(np.min([th1, th2]) - 1, np.max([th1, th2]) + 1)
ax2.set_ylim(np.min([w1, w2]) - 1, np.max([w1, w2]) + 1)
ax2.grid()

(traj1,) = ax2.plot(
    [], [], "-", lw=1, color="blue", alpha=0.6, label=r"$\theta_1 - \omega_1$"
)
(traj2,) = ax2.plot(
    [], [], "-", lw=1, color="red", alpha=0.6, label=r"$\theta_2 - \omega_2$"
)
(pt1,) = ax2.plot([], [], "o", color="blue")
(pt2,) = ax2.plot([], [], "o", color="red")
ax2.legend(loc="upper right")

history_len = 200  # 軌跡をどれくらい残すか


def update(i):
    # 左：実空間
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]
    line.set_data(thisx, thisy)
    mass1.set_data([x1[i]], [y1[i]])
    mass2.set_data([x2[i]], [y2[i]])

    start_idx = max(0, i - history_len)
    trace.set_data(x2[start_idx:i], y2[start_idx:i])

    time_text.set_text(f"Time = {t[i]:.2f} s")

    traj1.set_data(th1[:i], w1[:i])
    traj2.set_data(th2[:i], w2[:i])

    pt1.set_data([th1[i]], [w1[i]])
    pt2.set_data([th2[i]], [w2[i]])

    return line, mass1, mass2, trace, time_text, traj1, traj2, pt1, pt2


ani = animation.FuncAnimation(fig, update, frames=len(t), interval=dt * 1000, blit=True)

plt.show()

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint

# --- 無次元パラメータ設定 ---
# 物理パラメータの「比率」のみが重要になります
l = 1.0  # 長さの比 (L2 / L1)
mu = 1.0  # 質の比 (m2 / m1)

# シミュレーション設定
# 時間は無次元時間 tau = t * sqrt(g/L1)
tau_max = 500.0  # 無次元シミュレーション時間
d_tau = 0.05  # 無次元時間刻み

# 初期状態: [theta1, omega1, theta2, omega2]
# omega_dim = omega_phys * sqrt(L1/g))
init_state = [0.01, 0.0, 0.001, 0.0]
# init_state = [np.pi / 2 + 0.01, 0.0, np.pi + 0.001, 0.0]


# --- 運動方程式  ---
def double_pendulum_derivs_dimensionless(state, tau, l, mu):
    th1, w1, th2, w2 = state

    # g = 1, L1 = 1, m1 = 1 として扱える
    # m2 -> mu, L2 -> l と置き換える

    delta = th2 - th1
    cos_delta = np.cos(delta)
    sin_delta = np.sin(delta)

    # 分母の計算
    # den1 = (m1 + m2)*L1 - m2*L1*cos^2(delta)
    #      -> (1 + mu) - mu * cos^2(delta)
    den1 = (1 + mu) - mu * cos_delta * cos_delta
    den2 = l * den1

    dydt = np.zeros_like(state)
    dydt[0] = w1
    dydt[2] = w2

    # 分子1 (theta1の加速度に対応)
    # g -> 1, L1 -> 1, m2 -> mu, L2 -> l
    num1 = (
        mu * w1**2 * sin_delta * cos_delta
        + mu * 1 * np.sin(th2) * cos_delta  # g=1
        + mu * l * w2**2 * sin_delta
        - (1 + mu) * 1 * np.sin(th1)  # g=1
    )
    dydt[1] = num1 / den1

    # 分子2 (theta2の加速度に対応)
    num2 = -mu * l * w2**2 * sin_delta * cos_delta + (1 + mu) * (
        1 * np.sin(th1) * cos_delta  # g=1
        - 1 * w1**2 * sin_delta  # L1=1
        - 1 * np.sin(th2)  # g=1
    )
    dydt[3] = num2 / den2

    return dydt


# --- Tangent Dynamics (MLE & SALI用) ---
def tangent_dynamics_dim(augmented_state, tau, l, mu):
    state = augmented_state[:4]
    v1 = augmented_state[4:8]
    v2 = augmented_state[8:12]

    # 1. Base dynamics
    dstate = double_pendulum_derivs_dimensionless(state, tau, l, mu)

    # 2. Finite Difference Approximation for Tangent Vectors
    eps = 1e-8

    # For vector v1
    state_p1 = state + eps * v1
    dstate_p1 = double_pendulum_derivs_dimensionless(state_p1, tau, l, mu)
    dv1 = (dstate_p1 - dstate) / eps

    # For vector v2
    state_p2 = state + eps * v2
    dstate_p2 = double_pendulum_derivs_dimensionless(state_p2, tau, l, mu)
    dv2 = (dstate_p2 - dstate) / eps

    return np.concatenate([dstate, dv1, dv2])


# --- Main Simulation Loop ---
print("Simulating (Dimensionless)...")

t_eval = np.arange(0, tau_max, d_tau)
n_steps = len(t_eval)

# Initialize vectors
v1_init = np.random.randn(4)
v1_init /= np.linalg.norm(v1_init)

v2_init = np.random.randn(4)
v2_init -= np.dot(v2_init, v1_init) * v1_init
v2_init /= np.linalg.norm(v2_init)

current_aug_state = np.concatenate([init_state, v1_init, v2_init])

state_history = [init_state]
mle_history = []
sali_history = []
time_points = []
cum_log_norm = 0.0

for i in range(n_steps - 1):
    t_curr = t_eval[i]
    t_next = t_eval[i + 1]

    sol = odeint(
        tangent_dynamics_dim, current_aug_state, [t_curr, t_next], args=(l, mu)
    )
    next_aug = sol[-1]

    state_next = next_aug[:4]
    v1_next = next_aug[4:8]
    v2_next = next_aug[8:12]

    # MLE
    norm_v1 = np.linalg.norm(v1_next)
    cum_log_norm += np.log(norm_v1)
    mle = cum_log_norm / t_next
    v1_next /= norm_v1

    # SALI
    norm_v2 = np.linalg.norm(v2_next)
    v2_next /= norm_v2
    diff_norm = np.linalg.norm(v1_next - v2_next)
    sum_norm = np.linalg.norm(v1_next + v2_next)
    sali = min(diff_norm, sum_norm)

    state_history.append(state_next)
    mle_history.append(mle)
    sali_history.append(sali)
    time_points.append(t_next)

    current_aug_state = np.concatenate([state_next, v1_next, v2_next])

state_history = np.array(state_history)
th1 = state_history[:, 0]
w1 = state_history[:, 1]
th2 = state_history[:, 2]
w2 = state_history[:, 3]

# 座標変換 (長さは L1=1 で正規化されているためそのまま使用)
x1 = np.sin(th1)
y1 = -np.cos(th1)
x2 = x1 + l * np.sin(th2)
y2 = y1 - l * np.cos(th2)

# --- Plotting ---
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
ax_real = axes[0, 0]
ax_phase = axes[0, 1]
ax_mle = axes[1, 0]
ax_sali = axes[1, 1]

# 1. Real Space (Normalized Lengths)
ax_real.set_title("Real Space (Normalized by $L_1$)")
ax_real.set_xlabel("$x / L_1$")
ax_real.set_ylabel("$y / L_1$")
# 範囲は l (L2/L1) に依存
limit = (1.0 + l) * 1.1
ax_real.set_xlim(-limit, limit)
ax_real.set_ylim(-limit, limit)
ax_real.set_aspect("equal")
ax_real.grid()
(line,) = ax_real.plot([], [], "o-", lw=2, color="black")
(trace,) = ax_real.plot([], [], "-", lw=0.5, color="gray", alpha=0.5)

# 2. Phase Space
ax_phase.set_title("Phase Space")
ax_phase.set_xlabel(r"Angle $\theta$")
ax_phase.set_ylabel(r"Dim-less Angular Velocity $\dot{\theta}$")
ax_phase.plot(th1, w1, lw=0.5, alpha=0.6, label=r"$\theta_1$")
ax_phase.plot(th2, w2, lw=0.5, alpha=0.6, label=r"$\theta_2$")
ax_phase.legend()
ax_phase.grid()

# 3. MLE
ax_mle.set_title("MLE (Dimensionless)")
ax_mle.set_xlabel(r"Dimensionless Time $\tau$")
ax_mle.set_ylabel("MLE")
ax_mle.plot(time_points, mle_history, color="green", lw=1.5)
ax_mle.grid()

# 4. SALI
ax_sali.set_title("SALI")
ax_sali.set_xlabel(r"Dimensionless Time $\tau$")
ax_sali.set_ylabel("SALI (log scale)")
ax_sali.set_yscale("log")
ax_sali.plot(time_points, sali_history, color="purple", lw=1.5)
ax_sali.grid()

plt.tight_layout()


# Animation Function
def update(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]
    line.set_data(thisx, thisy)

    history_len = 200
    start_idx = max(0, i - history_len)
    trace.set_data(x2[start_idx:i], y2[start_idx:i])
    return line, trace


ani = animation.FuncAnimation(fig, update, frames=len(t_eval), interval=30, blit=True)
plt.show()

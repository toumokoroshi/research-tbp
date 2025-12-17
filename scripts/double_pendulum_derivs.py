import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint

# --- 無次元パラメータ設定 ---
l = 1.0  # 長さの比 (L2 / L1)
mu = 1.0  # 質量の比 (m2 / m1)

# シミュレーション設定
tau_max = 200.0  # カオスを見るため少し長めに
d_tau = 0.05
t_eval = np.arange(0, tau_max, d_tau)

# init_state = [0.0, 0.0, 0.1, 0.0]
init_state = [np.pi / 2 + 0.0001, 0, np.pi / 4, 0.0]

# 次元数
N_DIM = 4


# --- 運動方程式 (Hamiltonian系に由来) ---
def double_pendulum_derivs_dimensionless(state, tau, l, mu):
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


# --- Tangent Dynamics ---
def tangent_dynamics(augmented_state, tau, l, mu, n_vectors):
    state = augmented_state[:N_DIM]
    vectors = [
        augmented_state[N_DIM + i * N_DIM : N_DIM + (i + 1) * N_DIM]
        for i in range(n_vectors)
    ]

    dstate = double_pendulum_derivs_dimensionless(state, tau, l, mu)

    # 有限差分で接ベクトルを計算 (線形化方程式の代用)
    eps = 1e-8
    d_vectors = []
    for v in vectors:
        state_p = state + eps * v
        dstate_p = double_pendulum_derivs_dimensionless(state_p, tau, l, mu)
        dv = (dstate_p - dstate) / eps
        d_vectors.append(dv)

    return np.concatenate([dstate] + d_vectors)


# --- Helper: GALI Calculation ---
def compute_gali_svd(vectors):
    """特異値分解を用いて平行六面体の体積(GALI)を計算"""
    V = np.column_stack(vectors)
    # 特異値を計算
    s = np.linalg.svd(V, compute_uv=False)
    # GALI = 特異値の積
    return np.prod(s)


# --- Main Simulation ---
print("Simulating Chaos...")

# 追跡する偏差ベクトルの数 (GALI4まで見るので4本)
N_VECTORS = 4

# ランダムな初期偏差ベクトル (正規化のみ)
vectors_curr = []
for i in range(N_VECTORS):
    v = np.random.randn(N_DIM)
    v /= np.linalg.norm(v)
    vectors_curr.append(v)

current_aug_state = np.concatenate([init_state] + vectors_curr)

# 履歴格納用
state_history = [init_state]
mle_history = []
sali_history = []
gali4_history = []
time_points = []

cum_log_norm = 0.0

for i in range(len(t_eval) - 1):
    t_curr = t_eval[i]
    t_next = t_eval[i + 1]

    # 1ステップ積分
    sol = odeint(
        tangent_dynamics, current_aug_state, [t_curr, t_next], args=(l, mu, N_VECTORS)
    )
    next_aug = sol[-1]

    state_next = next_aug[:N_DIM]
    vectors_next = [
        next_aug[N_DIM + j * N_DIM : N_DIM + (j + 1) * N_DIM] for j in range(N_VECTORS)
    ]

    # --- ここが修正点 ---
    # 直交化(Gram-Schmidt)はせず、正規化(Normalize)だけ行う
    # ※ 直交化するとベクトルが強制的に開いてしまい、SALI/GALIが下がらなくなる

    vectors_normalized = []
    for j, v in enumerate(vectors_next):
        norm_v = np.linalg.norm(v)

        # MLEの計算 (1本目のベクトルの伸び率)
        if j == 0:
            cum_log_norm += np.log(norm_v)
            mle = cum_log_norm / t_next

        # オーバーフロー防止のための正規化 (向きは変えない！)
        if norm_v > 0:
            v_normed = v / norm_v
        else:
            v_normed = v
        vectors_normalized.append(v_normed)

    # SALI 計算 (v1 と v2 の平行度)
    v1 = vectors_normalized[0]
    v2 = vectors_normalized[1]
    diff_norm = np.linalg.norm(v1 - v2)
    sum_norm = np.linalg.norm(v1 + v2)
    sali = min(diff_norm, sum_norm)

    # GALI4 計算 (4本の相関)
    gali4 = compute_gali_svd(vectors_normalized[:4])

    # 履歴保存
    state_history.append(state_next)
    mle_history.append(mle)
    sali_history.append(sali)
    gali4_history.append(gali4)
    time_points.append(t_next)

    # 次のステップへ
    current_aug_state = np.concatenate([state_next] + vectors_normalized)

# --- Plotting ---
state_history = np.array(state_history)
th1, w1, th2, w2 = state_history.T

# 座標変換
x1 = np.sin(th1)
y1 = -np.cos(th1)
x2 = x1 + l * np.sin(th2)
y2 = y1 - l * np.cos(th2)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 1. Real Space Trajectory
ax_real = axes[0, 0]
ax_real.set_title(f"Real Space Trajectory (L={l}, mu={mu})")
ax_real.plot(x2, y2, lw=0.5, color="black", alpha=0.6)
ax_real.set_aspect("equal")
ax_real.grid()

# 2. MLE
ax_mle = axes[0, 1]
ax_mle.set_title("MLE (Lyapunov Exponent)")
ax_mle.plot(time_points, mle_history, color="green")
ax_mle.set_xlabel("Time")
ax_mle.set_ylabel("MLE")
ax_mle.grid()

# 3. SALI (Log Scale)
ax_sali = axes[1, 0]
ax_sali.set_title("SALI (Log Scale)")
ax_sali.semilogy(time_points, sali_history, color="purple")
ax_sali.set_xlabel("Time")
ax_sali.set_ylabel("SALI")
ax_sali.grid()

# 4. GALI4 (Log Scale)
ax_gali = axes[1, 1]
ax_gali.set_title("GALI 4 (Log Scale)")
ax_gali.semilogy(time_points, gali4_history, color="red")
ax_gali.set_xlabel("Time")
ax_gali.set_ylabel("GALI 4")
ax_gali.grid()

plt.tight_layout()
plt.show()

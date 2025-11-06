import numpy as np
from poliastro.bodies import Sun, Earth, Moon
# G (重力定数) の正確な値もインポート
from astropy.constants import G

# --- 1. 定数と質量の準備 ---
# poliastroからGM値 (k) を取得 (単位: km^3/s^2)
GM_sun = Sun.k.to_value('km**3/s**2')
GM_earth = Earth.k.to_value('km**3/s**2')
GM_moon = Moon.k.to_value('km**3/s**2')

# G (km^3 / (kg s^2))
G_val = G.to_value('km**3/(kg*s**2)')

# 質量を計算
m1 = GM_sun / G_val  # Sun
m2 = GM_earth / G_val # Earth
m3 = GM_moon / G_val # Moon

# --- 2. HORIZONSからのデータ入力 ---
# (例として2025-11-19 20:50 UT のデータを入力)
# !! 重要: ここの値は、必ずご自身でHORIZONSから取得した
# !! SSB基準, J2000 の値に置き換えてください !!

# SSB基準のJ2000慣性座標系 (単位: km, km/s)
R_sun = np.array([-6.65345719E+07, 1.29410313E+08, 5.61053073E+07])
V_sun = np.array([-2.48356872E+01, -1.33081533E+01, -5.76677840E+00])

R_earth = np.array([8.22152200E+07, 1.23112835E+08, 5.33718915E+07])
V_earth = np.array([-2.45788874E+01, 1.35010620E+01, 5.85246726E+00])

R_moon = np.array([7.865832465196744E+07, 1.133601172465368E+08, 4.912582040442912E+07])
V_moon = np.array([-2.482216875862968E+01, 1.407179045703669E+01, 6.066799863831345E+00])


# --- 3. CR3BPパラメータの計算 ---

# 質量パラメータ (μ)
mu = m2 / (m1 + m2)
print(f"CR3BP 質量パラメータ (μ = m2 / (m1 + m2)): {mu:.10e}")

# --- 4. 単位系の正規化 ---
# 太陽-地球の相対ベクトル (慣性系)
R_12 = R_earth - R_sun  # m1 -> m2 のベクトル
V_12 = V_earth - V_sun

# 単位長さ (L): 太陽と地球の距離
L = np.linalg.norm(R_12)
print(f"単位長さ (L): {L:.5e} km")

# システムの角運動量ベクトル (z軸の基準)
H_12 = np.cross(R_12, V_12)

# 単位時間 (T): 太陽-地球の公転周期に基づく
# n = sqrt(G(m1+m2) / L^3)
# T = 1 / n
n = np.sqrt( (GM_sun + GM_earth) / L**3 )
T = 1.0 / n
print(f"単位時間 (T): {T:.5e} s  (1/n)")

# 単位速度 (V = L/T)
V_unit = L / T
print(f"単位速度 (V = L/T): {V_unit:.5e} km/s")

# --- 5. 回転座標系の基底ベクトル (i, j, k) の定義 ---

# x軸 (i): 太陽から地球への方向
i_hat = R_12 / L

# z軸 (k): 系の角運動量の方向
k_hat = H_12 / np.linalg.norm(H_12)

# y軸 (j): 右手系 (j = k x i)
j_hat = np.cross(k_hat, i_hat)

# 回転行列 (慣性系 -> 回転系)
# Q[0,:] = i_hat, Q[1,:] = j_hat, Q[2,:] = k_hat
Q = np.array([i_hat, j_hat, k_hat])

# --- 6. 座標系の原点 (太陽-地球 重心) の計算 ---
R_bary = (m1 * R_sun + m2 * R_earth) / (m1 + m2)
V_bary = (m1 * V_sun + m2 * V_earth) / (m1 + m2)

# --- 7. 月のベクトルを変換 ---

# 1. 月のベクトルを「慣性系のまま、原点だけ重心に」移動
r_moon_inertial = R_moon - R_bary
v_moon_inertial = V_moon - V_bary

# 2. 回転座標系の角速度ベクトル (ω = n * k_hat)
# (CR3BPモデルでは n*k_hat を使う。現実の瞬時値 H_12/L**2 でも良い)
omega_vec = n * k_hat

# 3. 速度の変換 (v_inertial = v_rotating + ω x r)
# よって v_rotating = v_inertial - ω x r
# (v_rotating はまだ慣性系で表現されている)
v_moon_rotating_inertial_basis = v_moon_inertial - np.cross(omega_vec, r_moon_inertial)

# 4. 回転行列 Q を使って、回転座標系の基底に投影
#    (同時に正規化)
r_norm = Q @ r_moon_inertial / L
v_norm = Q @ v_moon_rotating_inertial_basis / V_unit


# --- 8. 結果の表示 ---
print("\n--- CR3BP シミュレータ用 初期値 ---")
print(f"正規化された月の位置 (x, y, z):")
print(f"  x = {r_norm[0]:.12f}")
print(f"  y = {r_norm[1]:.12f}")
print(f"  z = {r_norm[2]:.12f}")
print(f"\n正規化された月の速度 (vx, vy, vz):")
print(f"  vx = {v_norm[0]:.12f}")
print(f"  vy = {v_norm[1]:.12f}")
print(f"  vz = {v_norm[2]:.12f}")

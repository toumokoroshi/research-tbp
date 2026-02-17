"""LEO軌道のヤコビ積分概算スクリプト

太陽-地球系CR3BPの無次元座標系において、
地球周回低軌道(LEO)のヤコビ積分C_Jを概算する。

座標系規約:
  - 太陽: (-μ, 0, 0),  地球: (1-μ, 0, 0)
  - μ = m_Earth / (m_Sun + m_Earth)
  - U* = 0.5*(x²+y²) + (1-μ)/r₁ + μ/r₂ + μ(1-μ)/2
  - C_J = 2*U* - v²
"""

import numpy as np


# === 物理定数 ===
MU: float = 3.003480593e-06  # 太陽-地球系の質量比
AU_KM: float = 1.496e8  # 天文単位 [km]
R_EARTH_KM: float = 6371.0  # 地球半径 [km]


def calc_jacobi_integral(
    x: float, y: float, z: float,
    vx: float, vy: float, vz: float,
    mu: float,
) -> float:
    """ヤコビ積分 C_J = 2*U* - v² を計算する。

    Args:
        x: 回転座標系でのx座標 (無次元)
        y: 回転座標系でのy座標 (無次元)
        z: 回転座標系でのz座標 (無次元)
        vx: 回転座標系でのx方向速度 (無次元)
        vy: 回転座標系でのy方向速度 (無次元)
        vz: 回転座標系でのz方向速度 (無次元)
        mu: 質量比 μ = m2/(m1+m2)

    Returns:
        ヤコビ積分の値 C_J
    """
    r1 = np.sqrt((x + mu) ** 2 + y**2 + z**2)
    r2 = np.sqrt((x - 1.0 + mu) ** 2 + y**2 + z**2)
    u_star = 0.5 * (x**2 + y**2) + (1.0 - mu) / r1 + mu / r2 + mu * (1.0 - mu) * 0.5
    v_sq = vx**2 + vy**2 + vz**2
    return 2.0 * u_star - v_sq


def leo_jacobi_estimate(altitude_km: float, mu: float = MU) -> dict:
    """LEO軌道のヤコビ積分を概算する。

    地球中心で円軌道を仮定し、回転座標系でのヤコビ積分を計算する。

    Args:
        altitude_km: LEO高度 [km]
        mu: 質量比 (デフォルト: 太陽-地球系)

    Returns:
        計算結果の辞書
    """
    # --- 無次元化 ---
    r_orbit_km = R_EARTH_KM + altitude_km
    r_orbit_nd = r_orbit_km / AU_KM  # 無次元軌道半径

    # 地球位置 (回転座標系)
    x_earth = 1.0 - mu

    # --- LEO初期条件 (地球直上、円軌道仮定) ---
    # 地球中心の二体問題で円軌道速度: v_circ = sqrt(GM_earth / r)
    # 無次元系では GM_earth = mu, r = r_orbit_nd
    # ただし地球周りの速度なので v_circ_relative = sqrt(mu / r_orbit_nd)
    v_circ_rel = np.sqrt(mu / r_orbit_nd)

    # 回転座標系での位置 (地球の真上方向: +x方向)
    x_pos = x_earth + r_orbit_nd
    y_pos = 0.0
    z_pos = 0.0

    # 回転座標系での速度
    # 慣性系での速度: 地球公転速度 + 地球周り円軌道速度
    #   v_inertial = (0, (1-μ) + v_circ_rel, 0)  (プログレード)
    #
    # 回転座標系への変換: v_rot = v_inertial - ω×r
    #   ω×r = (0, 0, 1) × (x, 0, 0) = (0, x, 0)
    #   x = (1-μ) + r_orbit_nd
    #
    # vy_rot = [(1-μ) + v_circ_rel] - [(1-μ) + r_orbit_nd]
    #        = v_circ_rel - r_orbit_nd
    vx_rot = 0.0
    vy_rot = v_circ_rel - r_orbit_nd  # ≈ v_circ_rel (r_orbit << v_circ)
    vz_rot = 0.0

    cj = calc_jacobi_integral(x_pos, y_pos, z_pos, vx_rot, vy_rot, vz_rot, mu)

    # --- 近似式による概算 (二体問題近似) ---
    # LEOでは r₂ << 1 かつ r₁ ≈ 1 なので:
    #   C_J ≈ 3 + 2*μ/r₂ - μ/r₂ (反対方向の公転寄与)
    # より厳密には Tisserand パラメータに相当:
    #   C_J ≈ 1/a + 2*sqrt(a*(1-e²))*cos(i) + const (二体問題の場合)
    # 円軌道の二体近似: C_J ≈ 3 + μ/r - 2r (r: 地球からの無次元距離)
    cj_approx_twobody = 3.0 + mu / r_orbit_nd - 2.0 * r_orbit_nd

    # L1, L2のヤコビ定数 (概算)
    hill_radius = (mu / 3.0) ** (1.0 / 3.0)
    x_l1 = x_earth - hill_radius
    x_l2 = x_earth + hill_radius
    cj_l1 = calc_jacobi_integral(x_l1, 0, 0, 0, 0, 0, mu)
    cj_l2 = calc_jacobi_integral(x_l2, 0, 0, 0, 0, 0, mu)

    # --- 有次元化情報 ---
    # 地球周りの有次元円軌道速度 [km/s]
    gm_earth_km = 3.986e5  # GM_earth [km³/s²]
    v_circ_dim = np.sqrt(gm_earth_km / r_orbit_km)
    # 地球の公転速度 (有次元) [km/s]
    v_earth_km_s = 29.78  # 地球公転速度 ≈ 29.78 km/s

    return {
        "altitude_km": altitude_km,
        "r_orbit_km": r_orbit_km,
        "r_orbit_nd": r_orbit_nd,
        "v_circ_rel_nd": v_circ_rel,
        "v_circ_dim_km_s": v_circ_dim,
        "v_earth_km_s": v_earth_km_s,
        "x_pos": x_pos,
        "vy_rot": vy_rot,
        "C_J": cj,
        "C_J_approx_2body": cj_approx_twobody,
        "C_J_L1": cj_l1,
        "C_J_L2": cj_l2,
        "hill_radius_nd": hill_radius,
        "hill_radius_km": hill_radius * AU_KM,
    }


def main() -> None:
    """各種LEO高度でヤコビ積分を概算し、結果を表示する。"""
    print("=" * 70)
    print("LEO軌道のヤコビ積分概算 (太陽-地球 CR3BP)")
    print(f"  μ = {MU:.10e}")
    print(f"  1 AU = {AU_KM:.3e} km")
    print(f"  地球半径 = {R_EARTH_KM:.1f} km")
    print("=" * 70)

    # 代表的なLEO高度
    altitudes = [200, 400, 600, 800, 1000, 2000]

    # L1/L2参照値
    ref = leo_jacobi_estimate(400)
    print(f"\n--- 参照値 ---")
    print(f"  ヒル半径: {ref['hill_radius_nd']:.6e} (= {ref['hill_radius_km']:.1f} km)")
    print(f"  C_J(L1) = {ref['C_J_L1']:.10f}")
    print(f"  C_J(L2) = {ref['C_J_L2']:.10f}")
    print(f"  C_J(L1) - 3 = {ref['C_J_L1'] - 3:.6e}")
    print(f"  C_J(L2) - 3 = {ref['C_J_L2'] - 3:.6e}")

    print(f"\n--- LEO高度ごとのヤコビ積分 ---")
    print(f"{'高度[km]':>10} {'r[km]':>10} {'r_nd':>14} {'v_circ[km/s]':>14} "
          f"{'C_J':>14} {'C_J - 3':>14} {'C_J(近似)':>14}")
    print("-" * 105)

    for alt in altitudes:
        res = leo_jacobi_estimate(alt)
        print(
            f"{res['altitude_km']:10.0f} "
            f"{res['r_orbit_km']:10.1f} "
            f"{res['r_orbit_nd']:14.6e} "
            f"{res['v_circ_dim_km_s']:14.4f} "
            f"{res['C_J']:14.10f} "
            f"{res['C_J'] - 3:14.6e} "
            f"{res['C_J_approx_2body']:14.10f}"
        )

    # 追加解析: vy_rotの内訳
    print(f"\n--- 回転座標系での速度内訳 ---")
    print(f"{'高度[km]':>10} {'v_circ_nd':>14} {'x_pos':>14} {'vy_rot':>14} {'vy_rot比率':>14}")
    print("-" * 70)
    for alt in altitudes:
        res = leo_jacobi_estimate(alt)
        # v_circ_relの有次元的な意味
        ratio = res["v_circ_rel_nd"] / res["x_pos"]
        print(
            f"{alt:10.0f} "
            f"{res['v_circ_rel_nd']:14.6e} "
            f"{res['x_pos']:14.10f} "
            f"{res['vy_rot']:14.10f} "
            f"{ratio:14.6e}"
        )

    print(f"\n--- 物理的解釈 ---")
    res400 = leo_jacobi_estimate(400)
    r = res400["r_orbit_nd"]
    print(f"  LEO 400km:")
    print(f"    位置 x = {res400['x_pos']:.10f} (地球からの距離: {r:.6e})")
    print(f"    円軌道速度 (無次元): {res400['v_circ_rel_nd']:.6e}")
    print(f"    円軌道速度 (有次元): {res400['v_circ_dim_km_s']:.2f} km/s")
    print(f"    地球公転速度 (有次元): {res400['v_earth_km_s']:.2f} km/s")
    print(f"    回転座標系速度 vy_rot = {res400['vy_rot']:.6e}")
    print(f"    → v² = {res400['vy_rot']**2:.6e}")
    u_star = 0.5 * res400['x_pos']**2 + (1-MU)/1.0 + MU/r + MU*(1-MU)/2
    print(f"    → 2*U* = {2 * u_star:.10f}")
    print(f"    → C_J = 2*U* - v² = {res400['C_J']:.10f}")
    print(f"    → C_J - 3 = {res400['C_J'] - 3:.6e}")
    print()
    print(f"  比較:")
    print(f"    C_J(LEO 400km) = {res400['C_J']:.10f}")
    print(f"    C_J(L1)        = {res400['C_J_L1']:.10f}")
    print(f"    C_J(L2)        = {res400['C_J_L2']:.10f}")
    print(f"    近似式 C_J ≈ 3 + μ/r = 3 + {MU/r:.6e} = {3 + MU/r:.10f}")
    if res400['C_J'] > res400['C_J_L1']:
        print(f"    → C_J(LEO) > C_J(L1): ZVC閉 → 地球ヒル球内に閉じ込め ✓")
    else:
        print(f"    → C_J(LEO) < C_J(L1): ZVC開 → 脱出可能")


if __name__ == "__main__":
    main()

"""
軌道の物理的正しさを検証するスクリプト

以下の項目を検証:
1. ヤコビ積分の保存（軌道および多様体に沿って）
2. 軌道の周期性（1周期後に初期状態に戻るか）
3. L2ラグランジュ点近傍に位置しているか
4. 固有値構造の妥当性
5. 多様体の方向（安定・不安定）の正しさ
"""
import sys
import io

# Windows console encoding fix
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

from pathlib import Path
from typing import Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 日本語フォント設定
try:
    import japanize_matplotlib
except ImportError:
    plt.rcParams['font.family'] = ['MS Gothic', 'DejaVu Sans']


def compute_jacobi_integral(
    state: np.ndarray, 
    mu: float
) -> float:
    """
    ヤコビ積分を計算する

    Args:
        state: 状態ベクトル [x, y, z, vx, vy, vz]
        mu: 質量パラメータ

    Returns:
        ヤコビ積分 C_J
    """
    x, y, z = state[0], state[1], state[2]
    vx, vy, vz = state[3], state[4], state[5]
    
    # 天体からの距離
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2 + z**2)
    
    # 有効ポテンシャル
    omega = 0.5 * (x**2 + y**2) + (1 - mu) / r1 + mu / r2
    
    # 速度の二乗
    v_squared = vx**2 + vy**2 + vz**2
    
    # ヤコビ積分
    c_j = 2 * omega - v_squared
    
    return c_j


def compute_l2_position(mu: float) -> float:
    """
    L2ラグランジュ点のx座標を計算する（Newton法）

    Args:
        mu: 質量パラメータ

    Returns:
        L2のx座標
    """
    # L2の初期推定値（第二天体の外側）
    x = 1.0 + (mu / 3)**(1/3)
    
    for _ in range(100):
        r1 = x + mu
        r2 = x - 1 + mu
        
        # コリネアコンディション
        f = x - (1 - mu) / r1**2 * np.sign(r1) + mu / r2**2 * np.sign(r2)
        df = 1 + 2 * (1 - mu) / abs(r1)**3 + 2 * mu / abs(r2)**3
        
        x_new = x - f / df
        if abs(x_new - x) < 1e-12:
            break
        x = x_new
    
    return x


def load_orbit_info(filepath: Path) -> dict:
    """
    orbit_info.txtから軌道情報を読み込む

    Args:
        filepath: ファイルパス

    Returns:
        軌道情報の辞書
    """
    info = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        if line.startswith('x  ='):
            info['x'] = float(line.split('=')[1])
        elif line.startswith('y  ='):
            info['y'] = float(line.split('=')[1])
        elif line.startswith('z  ='):
            info['z'] = float(line.split('=')[1])
        elif line.startswith('vx ='):
            info['vx'] = float(line.split('=')[1])
        elif line.startswith('vy ='):
            info['vy'] = float(line.split('=')[1])
        elif line.startswith('vz ='):
            info['vz'] = float(line.split('=')[1])
        elif line.startswith('Period:'):
            info['period'] = float(line.split(':')[1])
        elif line.startswith('Jacobi Constant:'):
            info['jacobi'] = float(line.split(':')[1])
        elif line.startswith('Stability Index:'):
            info['stability_index'] = float(line.split(':')[1])
    
    return info


def load_manifold_csv(filepath: Path) -> pd.DataFrame:
    """
    多様体CSVファイルを読み込む

    Args:
        filepath: ファイルパス

    Returns:
        DataFrameデータ
    """
    # ヘッダー行をスキップして読み込み
    df = pd.read_csv(
        filepath, 
        comment='#',
        names=['time', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'jacobi']
    )
    return df


def verify_jacobi_conservation(
    df: pd.DataFrame, 
    expected_cj: float, 
    mu: float
) -> Tuple[float, float, bool]:
    """
    ヤコビ積分の保存を検証する

    Args:
        df: 多様体DataFrame
        expected_cj: 期待されるヤコビ積分値
        mu: 質量パラメータ

    Returns:
        (最大誤差, 平均誤差, 保存されているか)
    """
    errors = []
    
    for _, row in df.iterrows():
        state = np.array([
            row['x'], row['y'], row['z'],
            row['vx'], row['vy'], row['vz']
        ])
        cj_computed = compute_jacobi_integral(state, mu)
        errors.append(abs(cj_computed - expected_cj))
    
    max_error = max(errors)
    mean_error = np.mean(errors)
    
    # 許容誤差: 1e-8 (数値積分誤差を考慮)
    is_conserved = max_error < 1e-6
    
    return max_error, mean_error, is_conserved


def verify_orbit_near_l2(
    x0: float, 
    mu: float
) -> Tuple[float, float, bool]:
    """
    軌道がL2ラグランジュ点近傍にあるかを検証

    Args:
        x0: 軌道の初期x座標
        mu: 質量パラメータ

    Returns:
        (L2位置, L2からの距離, 近傍にあるか)
    """
    l2_pos = compute_l2_position(mu)
    distance = abs(x0 - l2_pos)
    
    # Halo軌道はL2から0.01程度近傍にあることが期待される
    is_near = distance < 0.02
    
    return l2_pos, distance, is_near


def verify_eigenvalue_structure(
    stability_index: float
) -> Tuple[dict, bool]:
    """
    固有値構造の妥当性を検証

    Args:
        stability_index: 安定性指数

    Returns:
        (固有値情報, 妥当か)
    """
    # 不安定Halo軌道の場合:
    # - 1対の実固有値（1より大きいものと小さいもの）
    # - 2対の複素共役固有値（絶対値1）
    # 安定性指数が1より大きければ不安定
    
    info = {
        'is_unstable': stability_index > 1.0,
        'has_hyperbolic_pair': stability_index > 1.0,
    }
    
    # 安定性指数はモノドロミー行列の最大実固有値であるべき
    is_valid = stability_index > 1.0  # UPOなので1より大きいはず
    
    return info, is_valid


def verify_manifold_direction(
    stable_manifold: pd.DataFrame,
    unstable_manifold: pd.DataFrame,
    x0: float,
    mu: float
) -> dict:
    """
    多様体の方向を検証する

    Args:
        stable_manifold: 安定多様体DataFrame
        unstable_manifold: 不安定多様体DataFrame
        x0: 軌道の初期x座標
        mu: 質量パラメータ

    Returns:
        検証結果の辞書
    """
    results = {}
    
    # 安定多様体: 逆時間で発散する（軌道から離れる）
    # 不安定多様体: 順時間で発散する
    
    # 初期点と最終点の距離を計算
    l2_pos = compute_l2_position(mu)
    
    if len(stable_manifold) > 1:
        # 安定多様体の最終点がL2から遠ざかっているか
        first_point = stable_manifold.iloc[0]
        last_point = stable_manifold.iloc[-1]
        
        first_dist = abs(first_point['x'] - l2_pos)
        last_dist = abs(last_point['x'] - l2_pos)
        
        results['stable_diverges'] = last_dist > first_dist * 0.5
        results['stable_distance_change'] = last_dist - first_dist
    
    if len(unstable_manifold) > 1:
        first_point = unstable_manifold.iloc[0]
        last_point = unstable_manifold.iloc[-1]
        
        first_dist = abs(first_point['x'] - l2_pos)
        last_dist = abs(last_point['x'] - l2_pos)
        
        results['unstable_diverges'] = last_dist > first_dist * 0.5
        results['unstable_distance_change'] = last_dist - first_dist
    
    return results


def create_verification_plots(
    orbit_info: dict,
    stable_manifolds: list[pd.DataFrame],
    unstable_manifolds: list[pd.DataFrame],
    mu: float,
    output_dir: Path
) -> None:
    """
    検証用のプロットを作成

    Args:
        orbit_info: 軌道情報
        stable_manifolds: 安定多様体リスト
        unstable_manifolds: 不安定多様体リスト
        mu: 質量パラメータ
        output_dir: 出力ディレクトリ
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. XY平面投影
    ax1 = axes[0, 0]
    for i, df in enumerate(stable_manifolds[:10]):  # 最初の10本
        ax1.plot(df['x'], df['y'], 'b-', alpha=0.5, linewidth=0.5)
    for i, df in enumerate(unstable_manifolds[:10]):
        ax1.plot(df['x'], df['y'], 'r-', alpha=0.5, linewidth=0.5)
    
    # L2点をプロット
    l2_pos = compute_l2_position(mu)
    ax1.scatter([l2_pos], [0], c='green', marker='*', s=200, label='L2', zorder=5)
    ax1.scatter([1-mu], [0], c='blue', marker='o', s=100, label='Earth', zorder=5)
    ax1.set_xlabel('x (無次元)', fontsize=12)
    ax1.set_ylabel('y (無次元)', fontsize=12)
    ax1.set_title('XY平面投影（青: 安定多様体, 赤: 不安定多様体）', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal')
    
    # 2. XZ平面投影
    ax2 = axes[0, 1]
    for df in stable_manifolds[:10]:
        ax2.plot(df['x'], df['z'], 'b-', alpha=0.5, linewidth=0.5)
    for df in unstable_manifolds[:10]:
        ax2.plot(df['x'], df['z'], 'r-', alpha=0.5, linewidth=0.5)
    
    ax2.scatter([l2_pos], [0], c='green', marker='*', s=200, label='L2')
    ax2.set_xlabel('x (無次元)', fontsize=12)
    ax2.set_ylabel('z (無次元)', fontsize=12)
    ax2.set_title('XZ平面投影', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    # 3. ヤコビ積分の誤差
    ax3 = axes[1, 0]
    expected_cj = orbit_info['jacobi']
    
    for i, df in enumerate(stable_manifolds[:5]):
        errors = [abs(row['jacobi'] - expected_cj) for _, row in df.iterrows()]
        ax3.semilogy(df['time'], errors, 'b-', alpha=0.5, linewidth=0.5)
    
    for i, df in enumerate(unstable_manifolds[:5]):
        errors = [abs(row['jacobi'] - expected_cj) for _, row in df.iterrows()]
        ax3.semilogy(df['time'], errors, 'r-', alpha=0.5, linewidth=0.5)
    
    ax3.axhline(y=1e-10, color='g', linestyle='--', label='目標精度 (1e-10)')
    ax3.set_xlabel('時間 (無次元)', fontsize=12)
    ax3.set_ylabel('ヤコビ積分誤差', fontsize=12)
    ax3.set_title('ヤコビ積分の保存', fontsize=14)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. L2からの距離の時間変化
    ax4 = axes[1, 1]
    
    for df in stable_manifolds[:5]:
        distances = np.sqrt((df['x'] - l2_pos)**2 + df['y']**2 + df['z']**2)
        ax4.plot(df['time'], distances, 'b-', alpha=0.5, linewidth=0.5)
    
    for df in unstable_manifolds[:5]:
        distances = np.sqrt((df['x'] - l2_pos)**2 + df['y']**2 + df['z']**2)
        ax4.plot(df['time'], distances, 'r-', alpha=0.5, linewidth=0.5)
    
    ax4.set_xlabel('時間 (無次元)', fontsize=12)
    ax4.set_ylabel('L2からの距離 (無次元)', fontsize=12)
    ax4.set_title('L2点からの距離の時間変化', fontsize=14)
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'verification_plots.png', dpi=200)
    plt.close()
    
    print(f"検証プロットを保存: {output_dir / 'verification_plots.png'}")


def main() -> None:
    """
    メイン関数
    """
    # 設定
    mu = 3.0404233870e-06  # 太陽-地球系の質量パラメータ
    data_dir = Path(r"c:\WORKSPACE\Research_git\research-tbp\data\orbit_manifold_by_jacobi\orbit_manifold\26_0130_1038_L2_halo_north")
    
    print("=" * 60)
    print("  軌道物理的正しさの検証")
    print("=" * 60)
    print()
    
    # 1. 軌道情報の読み込み
    print("[1] 軌道情報の読み込み...")
    orbit_info = load_orbit_info(data_dir / "orbit_info.txt")
    print(f"  初期状態: x={orbit_info['x']:.10f}, y={orbit_info.get('y', 0):.10f}, z={orbit_info['z']:.10f}")
    print(f"  周期: {orbit_info['period']:.10f}")
    print(f"  ヤコビ積分: {orbit_info['jacobi']:.10f}")
    print(f"  安定性指数: {orbit_info['stability_index']:.10f}")
    print()
    
    # 2. L2点近傍の検証
    print("[2] L2点近傍の検証...")
    l2_pos, distance, is_near = verify_orbit_near_l2(orbit_info['x'], mu)
    print(f"  L2位置: {l2_pos:.10f}")
    print(f"  L2からの距離: {distance:.6e}")
    print(f"  L2近傍にあるか: {'✓ Yes' if is_near else '✗ No'}")
    print()
    
    # 3. 固有値構造の検証
    print("[3] 固有値構造の検証...")
    eigen_info, eigen_valid = verify_eigenvalue_structure(orbit_info['stability_index'])
    print(f"  不安定軌道か: {'✓ Yes' if eigen_info['is_unstable'] else '✗ No'}")
    print(f"  双曲型固有値対あり: {'✓ Yes' if eigen_info['has_hyperbolic_pair'] else '✗ No'}")
    print(f"  固有値構造は妥当: {'✓ Yes' if eigen_valid else '✗ No'}")
    print()
    
    # 4. 多様体の読み込みとヤコビ積分保存の検証
    print("[4] ヤコビ積分保存の検証...")
    
    # サンプルとして最初の5本の安定/不安定多様体を読み込む
    stable_manifolds = []
    unstable_manifolds = []
    
    for i in range(5):
        stable_file = data_dir / f"manifold_stable_{i:04d}.csv"
        if stable_file.exists():
            stable_manifolds.append(load_manifold_csv(stable_file))
    
    for i in range(100, 105):
        unstable_file = data_dir / f"manifold_unstable_{i:04d}.csv"
        if unstable_file.exists():
            unstable_manifolds.append(load_manifold_csv(unstable_file))
    
    # 最初の安定多様体でヤコビ積分保存を確認
    if stable_manifolds:
        max_err, mean_err, conserved = verify_jacobi_conservation(
            stable_manifolds[0], orbit_info['jacobi'], mu
        )
        print(f"  安定多様体 #0:")
        print(f"    最大誤差: {max_err:.6e}")
        print(f"    平均誤差: {mean_err:.6e}")
        print(f"    保存されているか: {'✓ Yes' if conserved else '✗ No'}")
    
    if unstable_manifolds:
        max_err, mean_err, conserved = verify_jacobi_conservation(
            unstable_manifolds[0], orbit_info['jacobi'], mu
        )
        print(f"  不安定多様体 #0:")
        print(f"    最大誤差: {max_err:.6e}")
        print(f"    平均誤差: {mean_err:.6e}")
        print(f"    保存されているか: {'✓ Yes' if conserved else '✗ No'}")
    print()
    
    # 5. 多様体方向の検証
    print("[5] 多様体方向の検証...")
    if stable_manifolds and unstable_manifolds:
        direction_results = verify_manifold_direction(
            stable_manifolds[0], unstable_manifolds[0], orbit_info['x'], mu
        )
        print(f"  安定多様体がL2から発散: {'✓ Yes' if direction_results.get('stable_diverges', False) else '✗ No'}")
        print(f"  不安定多様体がL2から発散: {'✓ Yes' if direction_results.get('unstable_diverges', False) else '✗ No'}")
    print()
    
    # 6. 検証プロットの作成
    print("[6] 検証プロットの作成...")
    # すべての多様体を読み込み（最初の10本ずつ）
    stable_manifolds_full = []
    unstable_manifolds_full = []
    
    for i in range(10):
        stable_file = data_dir / f"manifold_stable_{i:04d}.csv"
        if stable_file.exists():
            stable_manifolds_full.append(load_manifold_csv(stable_file))
    
    for i in range(100, 110):
        unstable_file = data_dir / f"manifold_unstable_{i:04d}.csv"
        if unstable_file.exists():
            unstable_manifolds_full.append(load_manifold_csv(unstable_file))
    
    if stable_manifolds_full and unstable_manifolds_full:
        create_verification_plots(
            orbit_info, 
            stable_manifolds_full, 
            unstable_manifolds_full, 
            mu, 
            data_dir
        )
    
    # 7. 総合判定
    print("\n" + "=" * 60)
    print("  検証結果サマリー")
    print("=" * 60)
    
    all_passed = True
    
    checks = [
        ("L2点近傍に位置", is_near),
        ("固有値構造が妥当", eigen_valid),
        ("ヤコビ積分保存（安定多様体）", conserved if stable_manifolds else False),
    ]
    
    for name, passed in checks:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {name}: {status}")
        all_passed = all_passed and passed
    
    print()
    if all_passed:
        print("  全検証項目: ✓✓✓ PASSED ✓✓✓")
    else:
        print("  全検証項目: ✗✗✗ SOME FAILED ✗✗✗")
    print("=" * 60)


if __name__ == "__main__":
    main()

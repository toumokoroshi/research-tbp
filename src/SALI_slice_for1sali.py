import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from tkinter import filedialog, Tk
import numpy as np

# ファイル選択ダイアログを開く（Tkinter ウィンドウを非表示）
root = Tk()
root.withdraw()
file_path = filedialog.askopenfilename(filetypes=[("Text Files", "*.txt")])

if not file_path:
    print("No file selected. Exiting.")
    exit()

# データの読み込み
try:
    # ヘッダー情報の読み込みとデータの開始行検出
    with open(file_path, 'r') as file:
        lines = file.readlines()
        header_info = []
        data_start_index = None
        for i, line in enumerate(lines):
            if line.strip().startswith("mesh_num"):
                data_start_index = i + 1
                break
            header_info.append(line.strip())
    
    if data_start_index is None:
        raise ValueError("Data section not found in the file.")
    
    # データ部分の読み込み
    data = pd.read_csv(file_path, sep="\s+", header=None, skiprows=data_start_index, names=[
        "mesh_num", "time", "x", "y", "z", "jacobi constant", "SALI"
    ])
except Exception as e:
    print(f"Error reading file: {e}")
    exit()

# ヘッダー情報の表示
print("Header Information:")
for line in header_info:
    print(line)

# データのプレビュー
print("\nData Preview:")
print(data.head())
print(data.info())

# 必要な列を数値型に変換
numeric_columns = ["x", "y", "z", "SALI"]
for col in numeric_columns:
    data[col] = pd.to_numeric(data[col], errors="coerce")

# 欠損値を含む行を削除
data = data.dropna(subset=numeric_columns)

if data.empty:
    print("Error: The dataset is empty or contains invalid data.")
    exit()

# 初期設定
z_min, z_max = data["z"].min(), data["z"].max()
z_initial = (z_min + z_max) / 2  # 初期値はzの中央値
tolerance = 0.0 # z値のスライス範囲

# SALIの全体的な最小値と最大値を計算
sali_min, sali_max = -1.41421, 1.41421  # -sqrt(2) から sqrt(2) までの範囲

def get_slice(z_target):
    """指定したz値のスライスを取得"""
    return data[(data["z"] >= z_target - tolerance) & (data["z"] <= z_target + tolerance)]

# プロットのセットアップ
fig, ax = plt.subplots(figsize=(10, 8))
fig.subplots_adjust(bottom=0.15)  # 下部余白を広げる
slice_data = get_slice(z_initial)
sc = ax.scatter(slice_data["x"], slice_data["y"], c=slice_data["SALI"], cmap="seismic", vmin=sali_min, vmax=sali_max, s=30)
cbar = plt.colorbar(sc, label="SALI")
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
title = ax.set_title(f"Z = {z_initial:.6f}")

# スライダーの追加
ax_slider = plt.axes([0.15, 0.05, 0.6, 0.03], facecolor="lightgoldenrodyellow") 
slider = Slider(ax_slider, "Z", z_min, z_max, valinit=z_initial, valstep=0.000001)

# スライダーの更新イベント
def update(val):
    """スライダーの値に応じてプロットを更新"""
    z_target = slider.val
    slice_data = get_slice(z_target)
    if not slice_data.empty:
        sc.set_offsets(slice_data[["x", "y"]].values)
        sc.set_array(slice_data["SALI"].values)
        # カラーバーの範囲は固定したままにする
        sc.set_clim(sali_min, sali_max)
    else:
        sc.set_offsets(np.empty((0, 2)))
        sc.set_array(np.array([]))
    title.set_text(f"Z = {z_target:.6f}")
    fig.canvas.draw_idle()

slider.on_changed(update)

# plt.tight_layout()
plt.show()

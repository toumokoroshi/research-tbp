from pathlib import Path
import re
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
import pandas as pd

FALLBACK_MU = 3.0034805945421924e-06
ASTRO_PARAM_PATH = Path(__file__).resolve().parents[1] / 'configs' / 'astro_param' / 'astro_param.txt'
ZERO_VELOCITY_COLOR = '#d35400'
FORBIDDEN_REGION_COLOR = '#ffe8d6'
JACOBI_PATTERN = re.compile(
    r'JACOBI(?:\s+(?:INTEGRAL|CONSTANT))?[^:=\n]*[:=]\s*([+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)',
    re.IGNORECASE,
)


def load_mu_parameter():
    gm_earth = None
    gm_sun = None
    try:
        with open(ASTRO_PARAM_PATH, 'r') as f:
            for raw_line in f:
                line = raw_line.split('//', 1)[0].strip()
                if not line or '=' not in line:
                    continue
                key, value = line.split('=', 1)
                try:
                    numeric_value = float(value.strip())
                except ValueError:
                    continue
                key_lower = key.strip().lower()
                if key_lower == 'gm_earth':
                    gm_earth = numeric_value
                elif key_lower == 'gm_sun':
                    gm_sun = numeric_value
        if gm_earth is not None and gm_sun is not None:
            return gm_earth / (gm_earth + gm_sun)
    except OSError:
        pass
    return FALLBACK_MU


MU_PARAMETER = load_mu_parameter()

try:
    # Windowsの場合
    plt.rcParams['font.family'] = 'Yu Gothic'
    # Macの場合 (どちらかの行を有効にする)
    # plt.rcParams['font.family'] = 'Hiragino Sans'
except:
    # 上記フォントがない場合のフォールバック
    plt.rcParams['font.family'] = 'sans-serif'

plt.rcParams['axes.unicode_minus'] = False
class SALIContourApp:
    def __init__(self, root):
        self.root = root
        self.root.title("SALIカラーコンター作成ツール")
        self.root.geometry("1200x800")

        self.data = None
        self.jacobi_constant = None
        self.mu = MU_PARAMETER
        self.create_widgets()

    def create_widgets(self):
        # 制御パネル
        control_frame = ttk.Frame(self.root, padding="10")
        control_frame.pack(side=tk.LEFT, fill=tk.Y)

        # ... (ファイル選択、ヘッダ行数、区切り文字、自動検出 の部分は変更なし) ...

        # ファイル選択
        ttk.Label(control_frame, text="ファイル選択:").pack(anchor=tk.W, pady=(0, 5))
        file_frame = ttk.Frame(control_frame)
        file_frame.pack(fill=tk.X, pady=(0, 10))

        self.file_path = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.file_path, width=30).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(file_frame, text="参照", command=self.browse_file).pack(side=tk.LEFT, padx=(5, 0))

        # ヘッダ行数
        ttk.Label(control_frame, text="ヘッダ行数:").pack(anchor=tk.W, pady=(0, 5))
        self.header_lines = tk.IntVar(value=0)
        ttk.Spinbox(control_frame, from_=0, to=100, textvariable=self.header_lines, width=10).pack(anchor=tk.W, pady=(0, 10))

        # 区切り文字
        ttk.Label(control_frame, text="区切り文字:").pack(anchor=tk.W, pady=(0, 5))
        self.delimiter = tk.StringVar(value="スペース")
        delimiter_combo = ttk.Combobox(control_frame, textvariable=self.delimiter,
                                       values=["スペース", "タブ", "コンマ", "セミコロン"], width=15)
        delimiter_combo.pack(anchor=tk.W, pady=(0, 10))

        # 自動ヘッダ検出
        self.auto_detect = tk.BooleanVar(value=True)
        ttk.Checkbutton(control_frame, text="ヘッダを自動検出",
                        variable=self.auto_detect).pack(anchor=tk.W, pady=(0, 10))


        hill_frame = ttk.LabelFrame(control_frame, text="??Hill???????")
        hill_frame.pack(fill=tk.X, pady=(0, 10))
        hill_frame.columnconfigure(1, weight=1)

        self.show_hill = tk.BooleanVar(value=False)
        ttk.Checkbutton(hill_frame, text="Hill????", variable=self.show_hill).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=2)

        ttk.Label(hill_frame, text="?? (AU)").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.hill_radius = tk.DoubleVar(value=0.01)
        ttk.Entry(hill_frame, textvariable=self.hill_radius, width=10).grid(row=1, column=1, sticky=tk.EW, pady=2)

        ttk.Label(hill_frame, text="??X (AU)").grid(row=2, column=0, sticky=tk.W, pady=2)
        self.hill_center_x = tk.DoubleVar(value=0.0)
        ttk.Entry(hill_frame, textvariable=self.hill_center_x, width=10).grid(row=2, column=1, sticky=tk.EW, pady=2)

        ttk.Label(hill_frame, text="??Y (AU)").grid(row=3, column=0, sticky=tk.W, pady=2)
        self.hill_center_y = tk.DoubleVar(value=0.0)
        ttk.Entry(hill_frame, textvariable=self.hill_center_y, width=10).grid(row=3, column=1, sticky=tk.EW, pady=2)
        # --- Z値スライス (スライダーに変更) ---
        ttk.Label(control_frame, text="Z値スライス (AU):").pack(anchor=tk.W, pady=(0, 5))

        self.z_value = tk.DoubleVar(value=0.0)
        self.z_tolerance = tk.DoubleVar(value=0.1)

        # スライダー用のフレーム
        slider_frame = ttk.Frame(control_frame)
        slider_frame.pack(fill=tk.X, pady=(0, 5))

        ttk.Label(slider_frame, text="Z =").pack(side=tk.LEFT)
        # スライダーの値表示ラベル
        self.z_label = ttk.Label(slider_frame, text=f"{self.z_value.get():.3f}", width=6)
        self.z_label.pack(side=tk.LEFT, padx=5)

        self.z_slider = ttk.Scale(slider_frame, from_=0.0, to=0.0, # 範囲はload_dataで設定
                                  variable=self.z_value,
                                  orient=tk.HORIZONTAL,
                                  command=self.update_z_label,
                                  state=tk.DISABLED) # 初期状態は無効
        self.z_slider.pack(side=tk.LEFT, fill=tk.X, expand=True)

        # 許容範囲 (±) 用のフレーム
        tolerance_frame = ttk.Frame(control_frame)
        tolerance_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(tolerance_frame, text="許容範囲 ±").pack(side=tk.LEFT)
        ttk.Entry(tolerance_frame, textvariable=self.z_tolerance, width=10).pack(side=tk.LEFT, padx=5)
        ttk.Label(tolerance_frame, text="AU").pack(side=tk.LEFT)
        # --- 変更ここまで ---

        # 実行ボタン
        ttk.Button(control_frame, text="データ読み込み",
                   command=self.load_data).pack(fill=tk.X, pady=(10, 5))
        ttk.Button(control_frame, text="カラーコンター作成",
                   command=self.create_contour).pack(fill=tk.X, pady=(0, 10))

        # データ情報表示
        ttk.Separator(control_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        ttk.Label(control_frame, text="データ情報:").pack(anchor=tk.W, pady=(0, 5))
        self.info_text = tk.Text(control_frame, height=10, width=35)
        self.info_text.pack(fill=tk.BOTH, expand=True)

        # プロットエリア
        plot_frame = ttk.Frame(self.root)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.figure, self.ax = plt.subplots(figsize=(8, 6))
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    def update_z_label(self, val):
        """スライダーの値が変更されたときにラベルを更新"""
        self.z_label.config(text=f"{float(val):.3f}")
    def update_slider_resolution(self, *args):
        """許容範囲(z_tolerance)の変更をスライダーの解像度(resolution)に反映"""
        try:
            # 許容範囲に入力された値を取得
            resolution = self.z_tolerance.get()

            if resolution <= 0:
                # 解像度が0または負の値になるとエラーになるため、
                # その場合は非常に小さな値（ほぼ連続）にフォールバック
                resolution = 1e-9

            if hasattr(self, 'z_slider'): # z_sliderが作成された後か確認
                self.z_slider.config(resolution=resolution)

        except (tk.TclError, ValueError):
            # Entryが一時的に空（""）になったり、不正な値（"abc"）になると
            # get() が失敗することがあるため、例外を無視します。
            pass
    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="データファイルを選択",
            filetypes=[("Text files", "*.txt *.dat *.csv"), ("All files", "*.*")]
        )
        if filename:
            self.file_path.set(filename)

    def get_delimiter(self):
        delimiter_map = {
            "スペース": None,  # Pandasのdelim_whitespaceを使用
            "タブ": "\t",
            "コンマ": ",",
            "セミコロン": ";"
        }
        return delimiter_map[self.delimiter.get()]

    def detect_header_lines(self, filepath, delimiter):
        """ヘッダ行数を自動検出"""
        with open(filepath, 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line:
                    continue
                # 数値データのみの行を探す
                if delimiter is None:
                    parts = line.split()
                else:
                    parts = line.split(delimiter)

                try:
                    # 最初の列が数値かチェック
                    float(parts[0])
                    return i  # この行からデータが始まる
                except (ValueError, IndexError):
                    continue
        return 0

    def extract_jacobi_constant(self, filepath, header_lines):
        """ヘッダーに記されたヤコビ積分を抽出"""
        lines_to_check = header_lines if header_lines and header_lines > 0 else 30
        try:
            with open(filepath, 'r') as f:
                for idx, line in enumerate(f):
                    if lines_to_check and idx >= lines_to_check:
                        break
                    match = JACOBI_PATTERN.search(line)
                    if match:
                        try:
                            return float(match.group(1))
                        except ValueError:
                            continue
        except OSError:
            return None
        return None

    def load_data(self):
        if not self.file_path.get():
            messagebox.showerror("エラー", "ファイルを選択してください")
            return

        self.jacobi_constant = None
        try:
            delimiter = self.get_delimiter()

            # ヘッダ行数の決定
            if self.auto_detect.get():
                header_lines = self.detect_header_lines(self.file_path.get(), delimiter)
                self.header_lines.set(header_lines)
            else:
                header_lines = self.header_lines.get()

            self.jacobi_constant = self.extract_jacobi_constant(self.file_path.get(), header_lines)

            # データ読み込み
            if delimiter is None:
                self.data = pd.read_csv(self.file_path.get(),
                                        skiprows=header_lines,
                                        delim_whitespace=True,
                                        header=None,
                                        names=['meshnum', 'SALI', 'x', 'y', 'z', 'vx', 'vy', 'vz'])
            else:
                self.data = pd.read_csv(self.file_path.get(),
                                        skiprows=header_lines,
                                        delimiter=delimiter,
                                        header=None,
                                        names=['meshnum', 'SALI', 'x', 'y', 'z', 'vx', 'vy', 'vz'])

            # --- スライダー設定の更新 ---
            z_min = self.data['z'].min()
            z_max = self.data['z'].max()

            # スライダーの範囲と状態を更新
            self.z_slider.config(from_=z_min, to=z_max, state=tk.NORMAL)
            self.update_slider_resolution()

            # スライダーの初期値をZ範囲の中央に設定
            default_z = (z_min + z_max) / 2
            self.z_value.set(default_z)
            self.update_z_label(default_z) # ラベルも更新
            # --- 更新ここまで ---

            # データ情報表示
            info = f"データ読み込み成功\n\n"
            info += f"総データ数: {len(self.data)}\n"
            info += f"ヘッダ行数: {header_lines}\n\n"
            info += f"SALI範囲: [{self.data['SALI'].min():.4f}, {self.data['SALI'].max():.4f}]\n"
            info += f"X範囲: [{self.data['x'].min():.4f}, {self.data['x'].max():.4f}] AU\n"
            info += f"Y範囲: [{self.data['y'].min():.4f}, {self.data['y'].max():.4f}] AU\n"
            info += f"Z範囲: [{self.data['z'].min():.4f}, {self.data['z'].max():.4f}] AU\n"
            if self.jacobi_constant is not None:
                info += f"Jacobi integral (C): {self.jacobi_constant:.6f}\n"
            else:
                info += "Jacobi integral (C): N/A\n"

            self.info_text.delete(1.0, tk.END)
            self.info_text.insert(1.0, info)

            messagebox.showinfo("成功", "データの読み込みが完了しました")

        except Exception as e:
            messagebox.showerror("エラー", f"データ読み込みエラー:\n{str(e)}")

    def draw_zero_velocity_curve(self, z_target, reference_data=None):
        """ヘッダーで取得したヤコビ積分からゼロ速度曲線を描画する"""
        if self.jacobi_constant is None or self.mu is None:
            return None

        data_ref = reference_data if reference_data is not None else self.data
        if data_ref is None or data_ref.empty:
            return None

        try:
            x_min = float(data_ref['x'].min())
            x_max = float(data_ref['x'].max())
            y_min = float(data_ref['y'].min())
            y_max = float(data_ref['y'].max())
        except (KeyError, TypeError, ValueError):
            return None

        x_span = max(x_max - x_min, 1e-3)
        y_span = max(y_max - y_min, 1e-3)
        margin = max(0.05, 0.15 * max(x_span, y_span))

        x_vals = np.linspace(x_min - margin, x_max + margin, 400)
        y_vals = np.linspace(y_min - margin, y_max + margin, 400)
        X, Y = np.meshgrid(x_vals, y_vals)

        z_plane = float(z_target)
        zvc = np.ma.masked_invalid(self._compute_zero_velocity_field(X, Y, z_plane))

        if isinstance(zvc, np.ma.MaskedArray):
            if zvc.count() == 0:
                return None
            min_val = float(zvc.min())
            max_val = float(zvc.max())
        else:
            finite_mask = np.isfinite(zvc)
            if not finite_mask.any():
                return None
            min_val = float(np.min(zvc[finite_mask]))
            max_val = float(np.max(zvc[finite_mask]))

        if min_val < 0:
            self.ax.contourf(
                X,
                Y,
                zvc,
                levels=[min_val, 0],
                colors=[FORBIDDEN_REGION_COLOR],
                alpha=0.35,
                zorder=0.5,
            )

        if min_val <= 0 <= max_val:
            contour = self.ax.contour(
                X,
                Y,
                zvc,
                levels=[0],
                colors=[ZERO_VELOCITY_COLOR],
                linewidths=2.0,
                zorder=3,
            )
            if contour.collections:
                return Line2D(
                    [],
                    [],
                    color=ZERO_VELOCITY_COLOR,
                    linewidth=2.0,
                    label=f"Zero-velocity curve (C={self.jacobi_constant:.4f})",
                )
        return None

    def _compute_zero_velocity_field(self, x_grid, y_grid, z_plane):
        """2Ω - C（ゼロ速度面）を計算"""
        mu = self.mu if self.mu is not None else MU_PARAMETER
        eps = 1e-9
        r1 = np.sqrt((x_grid + mu) ** 2 + y_grid ** 2 + z_plane ** 2)
        r2 = np.sqrt((x_grid - (1 - mu)) ** 2 + y_grid ** 2 + z_plane ** 2)
        r1 = np.maximum(r1, eps)
        r2 = np.maximum(r2, eps)
        omega = 0.5 * (x_grid ** 2 + y_grid ** 2) + (1 - mu) / r1 + mu / r2
        return 2 * omega - self.jacobi_constant

    def create_contour(self):
        if self.data is None:
            messagebox.showerror("エラー", "先にデータを読み込んでください")
            return

        try:
            # Z値でスライス
            z_target = self.z_value.get()
            z_tol = self.z_tolerance.get()

            sliced_data = self.data[
                (self.data['z'] >= z_target - z_tol) &
                (self.data['z'] <= z_target + z_tol)
            ].copy()

            if len(sliced_data) == 0:
                messagebox.showerror("エラー", f"Z = {z_target} ± {z_tol} の範囲にデータがありません")
                return

            # カスタムカラーマップの作成
            sqrt2 = np.sqrt(2)
            colors = ['blue', 'white', 'red']
            n_bins = 256
            positions = [0, 1 / (sqrt2 + 1), 1]
            cmap = LinearSegmentedColormap.from_list('sali_cmap',
                                                     list(zip(positions, colors)),
                                                     N=n_bins)

            # --- ▼▼▼ ここから修正 ▼▼▼ ---

            # 1. self.ax.clear() の代わりに self.figure.clear() を使う
            # self.ax.clear() # <- 削除
            self.figure.clear() # <- これに変更

            # 2. figureをクリアしたので、ax を再作成する
            self.ax = self.figure.add_subplot(1, 1, 1)

            # --- ▲▲▲ 修正ここまで ---

            # プロット
            # 散布図として表示（補間なし）
            scatter = self.ax.scatter(sliced_data['x'], sliced_data['y'],
                                      c=sliced_data['SALI'],
                                      cmap=cmap,
                                      vmin=-1, vmax=sqrt2,
                                      s=20, marker='s', edgecolors='none', zorder=2)

            self.ax.set_xlabel('X (AU)', fontsize=12)
            self.ax.set_ylabel('Y (AU)', fontsize=12)
            self.ax.set_title(f'SALI カラーコンター (Z = {z_target:.3f} ± {z_tol:.3f} AU)',
                             fontsize=14)
            self.ax.set_aspect('equal')
            self.ax.grid(True, alpha=0.3)

            legend_handles = []
            if self.show_hill.get():
                try:
                    hill_radius = float(self.hill_radius.get())
                    hill_cx = float(self.hill_center_x.get())
                    hill_cy = float(self.hill_center_y.get())
                except (tk.TclError, ValueError):
                    messagebox.showerror("?????", "Hill???????????????????")
                    return

                if hill_radius <= 0:
                    messagebox.showerror("?????", "Hill??????????????????")
                    return

                hill_circle = Circle((hill_cx, hill_cy), hill_radius, fill=False, linestyle='--',
                                     linewidth=1.5, edgecolor='green', label='Earth Hill sphere')
                self.ax.add_patch(hill_circle)
                legend_handles.append(hill_circle)

            if legend_handles:
                self.ax.legend(handles=legend_handles, loc='upper right')

            # --- ▼▼▼ ここから修正 ▼▼▼ ---

            # 3. 古いカラーバーを削除する if ブロック全体を削除
            # if hasattr(self, 'colorbar'):
            #     self.colorbar.remove()

            # 4. 常に新しいカラーバーを作成する
            self.colorbar = self.figure.colorbar(scatter, ax=self.ax)

            # --- ▲▲▲ 修正ここまで ---

            self.colorbar.set_label('SALI', fontsize=12)

            self.figure.tight_layout()
            self.canvas.draw()

            # 情報更新
            # ... (変更なし) ...

        except Exception as e:
            messagebox.showerror("エラー", f"カラーコンター作成エラー:\n{str(e)}")

def main():
    root = tk.Tk()
    app = SALIContourApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()

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
ASTRO_PARAM_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "astro_param" / "astro_param.txt"
)
ZERO_VELOCITY_COLOR = "#d35400"
FORBIDDEN_REGION_COLOR = "#ffe8d6"
JACOBI_PATTERN = re.compile(
    r"JACOBI(?:\s+(?:INTEGRAL|CONSTANT))?[^:=\n]*[:=]\s*([+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)",
    re.IGNORECASE,
)


def load_mu_parameter():
    gm_earth = None
    gm_sun = None
    try:
        with open(ASTRO_PARAM_PATH, "r") as f:
            for raw_line in f:
                line = raw_line.split("//", 1)[0].strip()
                if not line or "=" not in line:
                    continue
                key, value = line.split("=", 1)
                try:
                    numeric_value = float(value.strip())
                except ValueError:
                    continue
                key_lower = key.strip().lower()
                if key_lower == "gm_earth":
                    gm_earth = numeric_value
                elif key_lower == "gm_sun":
                    gm_sun = numeric_value
        if gm_earth is not None and gm_sun is not None:
            return gm_earth / (gm_earth + gm_sun)
    except OSError:
        pass
    return FALLBACK_MU


MU_PARAMETER = load_mu_parameter()
print(f"Loaded mu parameter: {MU_PARAMETER:.10f}")

try:
    # Windowsの場合
    plt.rcParams["font.family"] = "Yu Gothic"
    # Macの場合 (どちらかの行を有効にする)
    # plt.rcParams['font.family'] = 'Hiragino Sans'
except:
    # 上記フォントがない場合のフォールバック
    plt.rcParams["font.family"] = "sans-serif"

plt.rcParams["axes.unicode_minus"] = False


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
        ttk.Entry(file_frame, textvariable=self.file_path, width=30).pack(
            side=tk.LEFT, fill=tk.X, expand=True
        )
        ttk.Button(file_frame, text="参照", command=self.browse_file).pack(
            side=tk.LEFT, padx=(5, 0)
        )

        # ヘッダ行数
        ttk.Label(control_frame, text="ヘッダ行数:").pack(anchor=tk.W, pady=(0, 5))
        self.header_lines = tk.IntVar(value=0)
        ttk.Spinbox(
            control_frame, from_=0, to=100, textvariable=self.header_lines, width=10
        ).pack(anchor=tk.W, pady=(0, 10))

        # 区切り文字
        ttk.Label(control_frame, text="区切り文字:").pack(anchor=tk.W, pady=(0, 5))
        self.delimiter = tk.StringVar(value="スペース")
        delimiter_combo = ttk.Combobox(
            control_frame,
            textvariable=self.delimiter,
            values=["スペース", "タブ", "コンマ", "セミコロン"],
            width=15,
        )
        delimiter_combo.pack(anchor=tk.W, pady=(0, 10))

        # 自動ヘッダ検出
        self.auto_detect = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            control_frame, text="ヘッダを自動検出", variable=self.auto_detect
        ).pack(anchor=tk.W, pady=(0, 10))

        hill_frame = ttk.LabelFrame(control_frame, text="??Hill???????")
        hill_frame.pack(fill=tk.X, pady=(0, 10))
        hill_frame.columnconfigure(1, weight=1)

        self.show_hill = tk.BooleanVar(value=False)
        ttk.Checkbutton(hill_frame, text="Hill????", variable=self.show_hill).grid(
            row=0, column=0, columnspan=2, sticky=tk.W, pady=2
        )

        ttk.Label(hill_frame, text="?? (AU)").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.hill_radius = tk.DoubleVar(value=0.01)
        ttk.Entry(hill_frame, textvariable=self.hill_radius, width=10).grid(
            row=1, column=1, sticky=tk.EW, pady=2
        )

        ttk.Label(hill_frame, text="??X (AU)").grid(
            row=2, column=0, sticky=tk.W, pady=2
        )

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
        ttk.Entry(file_frame, textvariable=self.file_path, width=30).pack(
            side=tk.LEFT, fill=tk.X, expand=True
        )
        ttk.Button(file_frame, text="参照", command=self.browse_file).pack(
            side=tk.LEFT, padx=(5, 0)
        )

        # ヘッダ行数
        ttk.Label(control_frame, text="ヘッダ行数:").pack(anchor=tk.W, pady=(0, 5))
        self.header_lines = tk.IntVar(value=0)
        ttk.Spinbox(
            control_frame, from_=0, to=100, textvariable=self.header_lines, width=10
        ).pack(anchor=tk.W, pady=(0, 10))

        # 区切り文字
        ttk.Label(control_frame, text="区切り文字:").pack(anchor=tk.W, pady=(0, 5))
        self.delimiter = tk.StringVar(value="スペース")
        delimiter_combo = ttk.Combobox(
            control_frame,
            textvariable=self.delimiter,
            values=["スペース", "タブ", "コンマ", "セミコロン"],
            width=15,
        )
        delimiter_combo.pack(anchor=tk.W, pady=(0, 10))

        # 自動ヘッダ検出
        self.auto_detect = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            control_frame, text="ヘッダを自動検出", variable=self.auto_detect
        ).pack(anchor=tk.W, pady=(0, 10))

        hill_frame = ttk.LabelFrame(control_frame, text="Hill球表示オプション")
        hill_frame.pack(fill=tk.X, pady=(0, 10))
        hill_frame.columnconfigure(1, weight=1)

        self.show_hill = tk.BooleanVar(value=False)
        ttk.Checkbutton(hill_frame, text="Hill球を表示", variable=self.show_hill).grid(
            row=0, column=0, columnspan=2, sticky=tk.W, pady=2
        )

        ttk.Label(hill_frame, text="半径 (AU)").grid(
            row=1, column=0, sticky=tk.W, pady=2
        )
        self.hill_radius = tk.DoubleVar(value=0.01)
        ttk.Entry(hill_frame, textvariable=self.hill_radius, width=10).grid(
            row=1, column=1, sticky=tk.EW, pady=2
        )

        ttk.Label(hill_frame, text="中心X (AU)").grid(
            row=2, column=0, sticky=tk.W, pady=2
        )
        self.hill_center_x = tk.DoubleVar(value=0.0)
        ttk.Entry(hill_frame, textvariable=self.hill_center_x, width=10).grid(
            row=2, column=1, sticky=tk.EW, pady=2
        )

        ttk.Label(hill_frame, text="中心Y (AU)").grid(
            row=3, column=0, sticky=tk.W, pady=2
        )
        self.hill_center_y = tk.DoubleVar(value=0.0)
        ttk.Entry(hill_frame, textvariable=self.hill_center_y, width=10).grid(
            row=3, column=1, sticky=tk.EW, pady=2
        )
        # --- 表示オプション ---
        ttk.Separator(hill_frame, orient=tk.HORIZONTAL).grid(
            row=4, column=0, columnspan=2, sticky=tk.EW, pady=5
        )

        self.show_collision = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            hill_frame, text="衝突(Collision)を表示", variable=self.show_collision
        ).grid(row=5, column=0, columnspan=2, sticky=tk.W, pady=2)

        self.show_escape = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            hill_frame, text="離脱(Escape)を表示", variable=self.show_escape
        ).grid(row=6, column=0, columnspan=2, sticky=tk.W, pady=2)

        # --- プロットモード選択 ---
        plot_mode_frame = ttk.LabelFrame(control_frame, text="プロットモード")
        plot_mode_frame.pack(fill=tk.X, pady=(0, 10))

        self.plot_mode = tk.StringVar(value="SALI")
        ttk.Radiobutton(
            plot_mode_frame, text="SALI", variable=self.plot_mode, value="SALI"
        ).pack(anchor=tk.W, pady=2)
        ttk.Radiobutton(
            plot_mode_frame,
            text="閾値到達時間 (lower_limit_reach_time)",
            variable=self.plot_mode,
            value="ReachTime",
        ).pack(anchor=tk.W, pady=2)

        # --- スケール設定 ---
        scale_frame = ttk.LabelFrame(control_frame, text="SALIスケール設定")
        scale_frame.pack(fill=tk.X, pady=(0, 10))

        self.scale_mode = tk.StringVar(value="Linear")
        ttk.Radiobutton(
            scale_frame, text="線形 (Linear)", variable=self.scale_mode, value="Linear"
        ).pack(anchor=tk.W, pady=2)
        ttk.Radiobutton(
            scale_frame, text="対数 (Log10)", variable=self.scale_mode, value="Log10"
        ).pack(anchor=tk.W, pady=2)

        # --- Z値スライス (Combobox選択式に変更) ---
        ttk.Label(control_frame, text="Z値スライス (AU):").pack(
            anchor=tk.W, pady=(0, 5)
        )

        # ユニークなZ値を保持する配列（load_dataで更新）
        self.z_values_unique = []
        self.z_value = tk.StringVar(value="")

        # Z値選択用のフレーム
        z_select_frame = ttk.Frame(control_frame)
        z_select_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(z_select_frame, text="Z =").pack(side=tk.LEFT)
        self.z_combo = ttk.Combobox(
            z_select_frame,
            textvariable=self.z_value,
            values=[],
            state="readonly",
            width=15,
        )
        self.z_combo.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)

        # Z値の数を表示するラベル
        self.z_count_label = ttk.Label(z_select_frame, text="(0 層)")
        self.z_count_label.pack(side=tk.LEFT)
        # --- 変更ここまで ---

        # 実行ボタン
        ttk.Button(control_frame, text="データ読み込み", command=self.load_data).pack(
            fill=tk.X, pady=(10, 5)
        )
        ttk.Button(
            control_frame, text="カラーコンター作成", command=self.create_contour
        ).pack(fill=tk.X, pady=(0, 10))

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

        # マウス座標表示用のラベル
        self.coord_label = ttk.Label(
            plot_frame, text="X: ---, Y: ---", font=("Consolas", 10)
        )
        self.coord_label.pack(anchor=tk.W, pady=(5, 0))

        # マウス移動イベントを接続
        self.canvas.mpl_connect("motion_notify_event", self.on_mouse_move)

    def on_mouse_move(self, event):
        """マウス移動時に座標を表示"""
        if event.inaxes == self.ax:
            self.coord_label.config(text=f"X: {event.xdata:.6f}, Y: {event.ydata:.6f}")
        else:
            self.coord_label.config(text="X: ---, Y: ---")

    def update_z_combo(self):
        """Z値のComboboxを更新"""
        if len(self.z_values_unique) > 0:
            # Z値を文字列リストに変換
            z_str_list = [f"{z:.15f}" for z in self.z_values_unique]
            self.z_combo.config(values=z_str_list)
            # デフォルトで最初の値を選択W
            self.z_value.set(z_str_list[0])
            self.z_count_label.config(text=f"({len(self.z_values_unique)} 層)")
        else:
            self.z_combo.config(values=[])
            self.z_value.set("")
            self.z_count_label.config(text="(0 層)")

    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="データファイルを選択",
            filetypes=[("Text files", "*.txt *.dat *.csv"), ("All files", "*.*")],
        )
        if filename:
            self.file_path.set(filename)

    def get_delimiter(self):
        delimiter_map = {
            "スペース": None,  # Pandasのdelim_whitespaceを使用
            "タブ": "\t",
            "コンマ": ",",
            "セミコロン": ";",
        }
        return delimiter_map[self.delimiter.get()]

    def detect_header_lines(self, filepath, delimiter):
        """ヘッダ行数を自動検出"""
        with open(filepath, "r") as f:
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

    def extract_header_info(self, filepath, header_lines):
        """ヘッダーから全ての情報を抽出"""
        header_info = {}
        if header_lines and header_lines > 0:
            lines_to_check = header_lines + 5
        else:
            lines_to_check = 30
        try:
            with open(filepath, "r") as f:
                raw_header_lines = []
                for idx, line in enumerate(f):
                    if lines_to_check and idx >= lines_to_check:
                        break
                    raw_header_lines.append(line.strip())

                    # Jacobi積分を抽出
                    match = JACOBI_PATTERN.search(line)
                    if match:
                        try:
                            header_info["jacobi_constant"] = float(match.group(1))
                        except ValueError:
                            pass

                    # KEY=VALUE または KEY: VALUE の形式を抽出
                    if "=" in line or ":" in line:
                        # コメント記号(#)の後を除去
                        clean_line = line.lstrip("#").strip()
                        for sep in ["=", ":"]:
                            if sep in clean_line:
                                parts = clean_line.split(sep, 1)
                                if len(parts) == 2:
                                    key = parts[0].strip().upper()
                                    value = parts[1].strip()
                                    # 既知のキーに対応
                                    if (
                                        "JACOBI" in key
                                        and "jacobi_constant" not in header_info
                                    ):
                                        try:
                                            header_info["jacobi_constant"] = float(
                                                value
                                            )
                                        except ValueError:
                                            pass
                                    elif "TIMESTEP" in key or "CALC_TIMESTEP" in key:
                                        header_info["timestep"] = value
                                    elif "TIME" in key and "THRESHOLD" in key:
                                        header_info["time_threshold"] = value
                                    elif "MESH" in key and (
                                        "NUM" in key or "COUNT" in key or "SIZE" in key
                                    ):
                                        header_info["mesh_info"] = value
                                    elif "MU" in key and key.replace(" ", "") in [
                                        "MU",
                                        "MUPARAMETER",
                                    ]:
                                        header_info["mu"] = value
                                    elif "CHAOS" in key or "INDEX" in key:
                                        header_info["chaos_index"] = value
                                    elif "SALI" in key and "LOWER" in key:
                                        header_info["sali_lower_limit"] = value
                                    elif "INCLINATION" in key:
                                        header_info["inclination"] = value
                                    elif "CONFIG" in key:
                                        header_info["config_file"] = value
                                break

                header_info["raw_header"] = raw_header_lines[
                    : min(15, len(raw_header_lines))
                ]
        except OSError:
            pass
        return header_info

    def extract_jacobi_constant(self, filepath, header_lines):
        """ヘッダーに記されたヤコビ積分を抽出（後方互換）"""
        info = self.extract_header_info(filepath, header_lines)
        return info.get("jacobi_constant")

    def load_data(self):
        if not self.file_path.get():
            messagebox.showerror("エラー", "ファイルを選択してください")
            return

        self.jacobi_constant = None
        self.header_info = {}
        try:
            delimiter = self.get_delimiter()

            # ヘッダ行数の決定
            if self.auto_detect.get():
                header_lines = self.detect_header_lines(self.file_path.get(), delimiter)
                self.header_lines.set(header_lines)
            else:
                header_lines = self.header_lines.get()

            # ヘッダ情報を抽出
            self.header_info = self.extract_header_info(
                self.file_path.get(), header_lines
            )
            self.jacobi_constant = self.header_info.get("jacobi_constant")

            # 列名の決定（自動判別）
            # まず1行だけ読んで列数を確認
            if delimiter is None:
                temp_df = pd.read_csv(
                    self.file_path.get(),
                    skiprows=header_lines,
                    sep=r"\s+",
                    header=None,
                    nrows=1,
                )
            else:
                temp_df = pd.read_csv(
                    self.file_path.get(),
                    skiprows=header_lines,
                    delimiter=delimiter,
                    header=None,
                    nrows=1,
                )

            num_cols = temp_df.shape[1]

            if num_cols == 11:
                col_names = [
                    "meshnum",
                    "SALI",
                    "x",
                    "y",
                    "z",
                    "vx",
                    "vy",
                    "vz",
                    "collision",
                    "escape",
                    "lower_limit_reach_time",
                ]
            elif num_cols == 10:
                col_names = [
                    "meshnum",
                    "SALI",
                    "x",
                    "y",
                    "z",
                    "vx",
                    "vy",
                    "vz",
                    "collision",
                    "escape",
                ]
            else:
                col_names = ["meshnum", "SALI", "x", "y", "z", "vx", "vy", "vz"]

            # データ読み込み
            if delimiter is None:
                self.data = pd.read_csv(
                    self.file_path.get(),
                    skiprows=header_lines,
                    sep=r"\s+",
                    header=None,
                    names=col_names,
                    low_memory=False,
                )
            else:
                self.data = pd.read_csv(
                    self.file_path.get(),
                    skiprows=header_lines,
                    delimiter=delimiter,
                    header=None,
                    names=col_names,
                    low_memory=False,
                )

            # データ型変換 (collision/escapeがない場合はNaNになるのでfillna等は不要だが、念のため)
            if "collision" not in self.data.columns:
                self.data["collision"] = 0
            if "escape" not in self.data.columns:
                self.data["escape"] = 0
            if "lower_limit_reach_time" not in self.data.columns:
                self.data["lower_limit_reach_time"] = -1  # 未到達を示すデフォルト値

            # --- Z値のユニーク値を抽出・ソート ---
            self.z_values_unique = np.sort(self.data["z"].unique())
            self.update_z_combo()
            # --- 更新ここまで ---

            # データ情報表示
            info = f"データ読み込み成功\n\n"
            info += f"=== ヘッダ情報 ===\n"
            if self.jacobi_constant is not None:
                info += f"Jacobi積分: {self.jacobi_constant:.6f}\n"
            if self.header_info.get("timestep"):
                info += f"タイムステップ: {self.header_info['timestep']}\n"
            if self.header_info.get("time_threshold"):
                info += f"積分時間: {self.header_info['time_threshold']}\n"
            if self.header_info.get("chaos_index"):
                info += f"カオス指標: {self.header_info['chaos_index']}\n"
            if self.header_info.get("sali_lower_limit"):
                info += f"SALI下限: {self.header_info['sali_lower_limit']}\n"
            if self.header_info.get("mu"):
                info += f"μ: {self.header_info['mu']}\n"
            if self.header_info.get("inclination"):
                info += f"傾斜角: {self.header_info['inclination']}\n"
            if self.header_info.get("mesh_info"):
                info += f"メッシュ: {self.header_info['mesh_info']}\n"
            if self.header_info.get("config_file"):
                info += f"設定ファイル: {self.header_info['config_file']}\n"

            info += f"\n=== データ統計 ===\n"
            info += f"総データ数: {len(self.data)}\n"
            info += f"ヘッダ行数: {header_lines}\n"
            info += f"SALI範囲: [{self.data['SALI'].min():.4f}, {self.data['SALI'].max():.4f}]\n"
            info += (
                f"X範囲: [{self.data['x'].min():.4f}, {self.data['x'].max():.4f}] AU\n"
            )
            info += (
                f"Y範囲: [{self.data['y'].min():.4f}, {self.data['y'].max():.4f}] AU\n"
            )
            info += (
                f"Z範囲: [{self.data['z'].min():.4f}, {self.data['z'].max():.4f}] AU\n"
            )

            # lower_limit_reach_time の範囲情報
            if "lower_limit_reach_time" in self.data.columns:
                valid_reach_time = self.data[self.data["lower_limit_reach_time"] >= 0][
                    "lower_limit_reach_time"
                ]
                if len(valid_reach_time) > 0:
                    info += f"\n閾値到達時間範囲: [{valid_reach_time.min():.4f}, {valid_reach_time.max():.4f}]\n"
                    info += f"(閾値到達データ数: {len(valid_reach_time)})\n"

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
            x_min = float(data_ref["x"].min())
            x_max = float(data_ref["x"].max())
            y_min = float(data_ref["y"].min())
            y_max = float(data_ref["y"].max())
        except (KeyError, TypeError, ValueError):
            return None

        x_span = max(x_max - x_min, 1e-3)
        y_span = max(y_max - y_min, 1e-3)
        margin = max(0.005, 0.05 * max(x_span, y_span))

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
            # Check if any contour lines were actually generated
            # Some versions of matplotlib might not populate collections immediately or behave differently
            has_contours = False
            if hasattr(contour, "collections"):
                if len(contour.collections) > 0:
                    has_contours = True
            elif hasattr(contour, "allsegs"):
                if len(contour.allsegs) > 0:
                    has_contours = True
            else:
                # Fallback: assume if we called contour, we want the legend
                has_contours = True

            if has_contours:
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
        r1 = np.sqrt((x_grid + mu) ** 2 + y_grid**2 + z_plane**2)
        r2 = np.sqrt((x_grid - (1 - mu)) ** 2 + y_grid**2 + z_plane**2)
        r1 = np.maximum(r1, eps)
        r2 = np.maximum(r2, eps)
        omega = 0.5 * (x_grid**2 + y_grid**2) + (1 - mu) / r1 + mu / r2
        return 2 * omega - self.jacobi_constant

    def create_contour(self):
        if self.data is None:
            messagebox.showerror("エラー", "先にデータを読み込んでください")
            return

        try:
            # Z値でスライス（選択したZ値に完全一致するデータを取得）
            z_str = self.z_value.get()
            if not z_str:
                messagebox.showerror("エラー", "Z値を選択してください")
                return

            z_target = float(z_str)

            # 浮動小数点精度の問題を回避するため np.isclose を使用
            sliced_data = self.data[
                np.isclose(self.data["z"], z_target, rtol=1e-9, atol=1e-12)
            ].copy()

            if len(sliced_data) == 0:
                messagebox.showerror("エラー", f"Z = {z_target} のデータがありません")
                return

            self.figure.clear()
            self.ax = self.figure.add_subplot(1, 1, 1)

            # プロットモードの取得
            plot_mode = self.plot_mode.get()

            legend_handles = []

            if plot_mode == "ReachTime":
                # --- 閾値到達時間モード ---
                # lower_limit_reach_time >= 0 のデータのみ有効
                valid_mask = sliced_data["lower_limit_reach_time"] >= 0
                valid_data = sliced_data[valid_mask]
                invalid_data = sliced_data[~valid_mask]

                # カラーマップ: viridis (青->緑->黄, 時間が長いほど明るい)
                cmap = "viridis"

                if len(valid_data) > 0:
                    plot_values = valid_data["lower_limit_reach_time"].values
                    vmin = plot_values.min()
                    vmax = plot_values.max()
                    label_str = "閾値到達時間 (無次元時間)"

                    scatter = self.ax.scatter(
                        valid_data["x"],
                        valid_data["y"],
                        c=plot_values,
                        cmap=cmap,
                        vmin=vmin,
                        vmax=vmax,
                        s=20,
                        marker="s",
                        edgecolors="none",
                        zorder=2,
                    )

                    # カラーバー
                    self.colorbar = self.figure.colorbar(scatter, ax=self.ax)
                    self.colorbar.set_label(label_str, fontsize=12)

                # 未到達データ（lower_limit_reach_time < 0）をグレーで表示
                if len(invalid_data) > 0:
                    self.ax.scatter(
                        invalid_data["x"],
                        invalid_data["y"],
                        c="lightgray",
                        s=20,
                        marker="s",
                        edgecolors="none",
                        zorder=1,
                        label="未到達 (閾値未達成)",
                    )

                title_str = f"閾値到達時間 カラーコンター (Z = {z_target:.6f} AU)"

            else:
                # --- SALIモード（デフォルト）---
                # SALI = -1 (計算不可) を分離
                invalid_sali_mask = sliced_data["SALI"] == -1
                valid_sali_mask = ~invalid_sali_mask

                valid_data = sliced_data[valid_sali_mask]
                invalid_data = sliced_data[invalid_sali_mask]

                # スケールモード
                scale_mode = self.scale_mode.get()

                # カラーマップ作成 (補色: 赤(0) <-> シアン(sqrt(2)))
                # 0(Chaos) -> Red, sqrt(2)(Order) -> Cyan
                colors = ["#ffffff", "#ff0000"]
                cmap = LinearSegmentedColormap.from_list(
                    "sali_complementary", colors, N=256
                )

                plot_values = valid_data["SALI"].values
                vmin = 0
                vmax = np.sqrt(2)

                if scale_mode == "Log10":
                    # Logスケールの場合、0以下の値は扱えないため、小さな値にクリップするかマスクする
                    # SALIは理論上0以上。0の場合は非常に小さな値として扱う
                    eps = 1e-16
                    plot_values = np.log10(np.maximum(plot_values, eps))
                    vmin = np.log10(eps)
                    vmax = np.log10(np.sqrt(2))
                    label_str = "Log10(SALI)"
                else:
                    label_str = "SALI"

                # --- プロット ---

                # 1. 無効なSALI (-1) -> 白 (背景と同化させるため、プロットしないか、白でプロット)
                if len(invalid_data) > 0:
                    self.ax.scatter(
                        invalid_data["x"],
                        invalid_data["y"],
                        c="blue",
                        s=20,
                        marker="s",
                        edgecolors="none",
                        zorder=1,
                        label="Undefined (SALI=-1)",
                    )

                # 2. 有効なSALI -> カラーマップ
                if len(valid_data) > 0:
                    scatter = self.ax.scatter(
                        valid_data["x"],
                        valid_data["y"],
                        c=plot_values,
                        cmap=cmap,
                        vmin=vmin,
                        vmax=vmax,
                        s=20,
                        marker="s",
                        edgecolors="none",
                        zorder=2,
                    )

                    # カラーバー
                    self.colorbar = self.figure.colorbar(scatter, ax=self.ax)
                    self.colorbar.set_label(label_str, fontsize=12)

                title_str = f"SALI カラーコンター (Z = {z_target:.6f} AU)"

            # 3. 衝突・離脱のオーバーレイ
            if self.show_collision.get():
                collision_data = sliced_data[sliced_data["collision"] == 1]
                if len(collision_data) > 0:
                    coll_scatter = self.ax.scatter(
                        collision_data["x"],
                        collision_data["y"],
                        c="black",
                        marker="x",
                        s=30,
                        linewidths=1,
                        zorder=4,
                        label="Collision",
                    )
                    legend_handles.append(coll_scatter)

            if self.show_escape.get():
                escape_data = sliced_data[sliced_data["escape"] == 1]
                if len(escape_data) > 0:
                    esc_scatter = self.ax.scatter(
                        escape_data["x"],
                        escape_data["y"],
                        c="black",
                        marker="^",
                        s=30,
                        linewidths=1,
                        zorder=4,
                        label="Escape",
                    )
                    legend_handles.append(esc_scatter)

            # --- 装飾 ---
            self.ax.set_xlabel("X (AU)", fontsize=12)
            self.ax.set_ylabel("Y (AU)", fontsize=12)
            self.ax.set_title(title_str, fontsize=14)
            self.ax.set_aspect("equal")
            self.ax.grid(True, alpha=0.3)

            if self.jacobi_constant is not None:
                self.ax.text(
                    0.02,
                    0.98,
                    f"C_J = {self.jacobi_constant:.4f}",
                    transform=self.ax.transAxes,
                    fontsize=10,
                    fontweight="bold",
                    verticalalignment="top",
                    bbox=dict(facecolor="white", alpha=0.5, edgecolor="none"),
                )

            zero_velocity_handle = self.draw_zero_velocity_curve(z_target)
            if zero_velocity_handle is not None:
                legend_handles.append(zero_velocity_handle)

            if self.show_hill.get():
                try:
                    hill_radius = float(self.hill_radius.get())
                    hill_cx = float(self.hill_center_x.get())
                    hill_cy = float(self.hill_center_y.get())
                except (tk.TclError, ValueError):
                    messagebox.showerror("エラー", "Hill球のパラメータが不正です")
                    return

                if hill_radius > 0:
                    hill_circle = Circle(
                        (hill_cx, hill_cy),
                        hill_radius,
                        fill=False,
                        linestyle="--",
                        linewidth=1.5,
                        edgecolor="green",
                        label="Earth Hill sphere",
                    )
                    self.ax.add_patch(hill_circle)
                    legend_handles.append(hill_circle)

            if legend_handles:
                self.ax.legend(handles=legend_handles, loc="upper right")

            self.figure.tight_layout()
            self.canvas.draw()

        except Exception as e:
            messagebox.showerror("エラー", f"カラーコンター作成エラー:\n{str(e)}")


def main():
    root = tk.Tk()
    app = SALIContourApp(root)

    def on_closing():
        """ウィンドウ閉じるときにアプリケーションを終了"""
        plt.close("all")  # matplotlibの全ウィンドウを閉じる
        root.destroy()
        root.quit()

    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.mainloop()


if __name__ == "__main__":
    main()

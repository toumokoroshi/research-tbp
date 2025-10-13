import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import customtkinter as ctk
from tkinter import filedialog, messagebox
import numpy as np
from scipy import stats

# 初期設定
ctk.set_appearance_mode("System")
ctk.set_default_color_theme("blue")

# グローバル変数
data = None
file_path = None
sali_min, sali_max = -1.41421, 1.41421
x_min, x_max, x_initial = None, None, None
unique_x_values = None


def load_file():
    global file_path, data, x_min, x_max, x_initial, unique_x_values

    file_path = filedialog.askopenfilename(
        filetypes=[("Text Files", "*.txt *.dat")])
    if not file_path:
        return

    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            data_start_index = next(i for i, line in enumerate(
                lines) if line.strip().startswith("mesh_num"))

        data = pd.read_csv(file_path, sep="\s+", header=None, skiprows=data_start_index, names=[
            "mesh_num", "time", "x", "y", "z", "jacobi constant", "SALI1", "SALI2", "SALI3", "SALI"
        ])

        numeric_columns = ["x", "y", "z",
                           "jacobi constant", "SALI1", "SALI2", "SALI3", "SALI"]
        for col in numeric_columns:
            data[col] = pd.to_numeric(data[col], errors="coerce")

        data.dropna(subset=numeric_columns, inplace=True)

        if data.empty:
            raise ValueError("The dataset is empty or contains invalid data.")
        # x値の初期化
        unique_x_values = sorted(data["x"].unique())
        x_min, x_max = min(unique_x_values), max(unique_x_values)
        x_initial = unique_x_values[len(unique_x_values) // 2]

        file_label.configure(text=f"Loaded File:\n{file_path}")
        messagebox.showinfo("Success", "File loaded successfully!")
    except Exception as e:
        messagebox.showerror("Error", f"Error reading file: {str(e)}")


def find_nearest_x(target):
    return unique_x_values[np.abs(np.array(unique_x_values) - target).argmin()]


def get_slice(x_target):
    return data[data["x"] == x_target]


def create_plot(column):
    if data is None:
        messagebox.showwarning(
            "Warning", "No data loaded. Please load a file first.")
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    fig.subplots_adjust(bottom=0.25, right=0.9)
    slice_data = get_slice(x_initial)

    # グリッドの作成
    y = slice_data["y"].values
    z = slice_data["z"].values
    sali = slice_data[column].values

    # グリッドポイントの作成
    yi = np.linspace(y.min(), y.max(), 100)
    zi = np.linspace(z.min(), z.max(), 100)
    yi, zi = np.meshgrid(yi, zi)

    # 補間
    from scipy.interpolate import griddata
    sali_i = griddata((y, z), sali, (yi, zi), method='cubic')

    # コンタープロット
    pc = ax.pcolormesh(yi, zi, sali_i, cmap=colormap_var.get(),
                       vmin=sali_min, vmax=sali_max, shading='auto')
    cbar = plt.colorbar(pc, label=column)

    ax.set_xlabel("Y-axis(au)")
    ax.set_ylabel("Z-axis(au)")
    title = ax.set_title(f"J = {slice_data['jacobi constant'].values[0]:.6f}")
    x_title = ax.text(
        0.5, 1.1, f"X = {x_initial:.6f}", transform=ax.transAxes, ha='center')

    ax_slider = plt.axes([0.2, 0.1, 0.6, 0.03],
                         facecolor="lightgoldenrodyellow")
    slider = Slider(ax_slider, "X", 0, len(unique_x_values) - 1,
                    valinit=unique_x_values.index(x_initial), valstep=1)

    def update(val):
        x_index = int(slider.val)
        x_target = unique_x_values[x_index]
        slice_data = get_slice(x_target)
        if not slice_data.empty:
            y = slice_data["y"].values
            z = slice_data["z"].values
            sali = slice_data[column].values
            sali_i = griddata((y, z), sali, (yi, zi), method='cubic')
            pc.set_array(sali_i.ravel())
            x_title.set_text(f"X = {x_target:.6f}")
        title.set_text(f"C_j = {slice_data['jacobi constant'].values[0]:.6f}")
        fig.canvas.draw_idle()

    slider.on_changed(update)

    # 座標表示の更新
    coord_display = ax.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points",
                                bbox=dict(boxstyle="round", fc="w"),
                                arrowprops=dict(arrowstyle="->"))
    coord_display.set_visible(False)

    def on_hover(event):
        if event.inaxes == ax:
            y, z = event.xdata, event.ydata
            coord_display.xy = (y, z)
            coord_display.set_text(f"({y:.8f}, {z:.8f})")
            coord_display.set_visible(True)
            fig.canvas.draw_idle()
        else:
            coord_display.set_visible(False)
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", on_hover)

    # ズーム機能は同じまま
    def on_zoom(event):
        if event.button == 'up':
            ax.set_xlim(ax.get_xlim()[0] * 0.9, ax.get_xlim()[1] * 0.9)
            ax.set_ylim(ax.get_ylim()[0] * 0.9, ax.get_ylim()[1] * 0.9)
        elif event.button == 'down':
            ax.set_xlim(ax.get_xlim()[0] * 1.1, ax.get_xlim()[1] * 1.1)
            ax.set_ylim(ax.get_ylim()[0] * 1.1, ax.get_ylim()[1] * 1.1)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect('scroll_event', on_zoom)

    plt.show()


def on_column_change(new_column):
    create_plot(new_column)


def show_help():
    help_text = """
    Usage Instructions:
    1. Click 'Select File' to load your data file.
    2. Choose the SALI column you want to visualize.
    3. Adjust marker size and color map if desired.
    4. Click 'Add Plot' to create the visualization.
    5. Use the slider to change the Z value.
    6. Hover over points to see exact coordinates.
    7. Use mouse wheel to zoom in/out.
    8. Click 'Save Plot' button to save the current plot as an image.
    """
    messagebox.showinfo("Help", help_text)


# GUI設定
app = ctk.CTk()
app.geometry("300x400")
app.title("SALI Plot Configuration")

file_label = ctk.CTkLabel(app, text="No file loaded",
                          wraplength=350, justify="left")
file_label.pack(pady=10)

file_button = ctk.CTkButton(app, text="Select File", command=load_file)
file_button.pack(pady=10)

label = ctk.CTkLabel(app, text="Select SALI Column:")
label.pack(pady=10)

column_var = ctk.StringVar(value="SALI1")
dropdown = ctk.CTkOptionMenu(app, variable=column_var, values=[
                             "SALI1", "SALI2", "SALI3", "SALI"], command=on_column_change)
dropdown.pack(pady=10)

# マーカーサイズ設定
marker_size_label = ctk.CTkLabel(app, text="Marker Size:")
marker_size_label.pack(pady=5)
marker_size_var = ctk.StringVar(value="50")
marker_size_entry = ctk.CTkEntry(app, textvariable=marker_size_var)
marker_size_entry.pack(pady=5)

# カラーマップ設定
colormap_label = ctk.CTkLabel(app, text="Color Map:")
colormap_label.pack(pady=5)
colormap_var = ctk.StringVar(value="seismic")
colormap_dropdown = ctk.CTkOptionMenu(app, variable=colormap_var, values=[
                                      "seismic", "viridis", "plasma", "inferno"])
colormap_dropdown.pack(pady=5)

add_plot_button = ctk.CTkButton(
    app, text="Add Plot", command=lambda: create_plot(column_var.get()))
add_plot_button.pack(pady=10)

stats_label = ctk.CTkLabel(app, text="", wraplength=350, justify="left")
stats_label.pack(pady=10)

help_button = ctk.CTkButton(app, text="Help", command=show_help)
help_button.pack(pady=10)

exit_button = ctk.CTkButton(app, text="Close", command=app.destroy)
exit_button.pack(pady=10)

app.mainloop()

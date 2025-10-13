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
z_min, z_max, z_initial = None, None, None
unique_z_values = None


def load_file():
    global file_path, data, z_min, z_max, z_initial, unique_z_values

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

        unique_z_values = sorted(data["z"].unique())
        z_min, z_max = min(unique_z_values), max(unique_z_values)
        z_initial = unique_z_values[len(unique_z_values) // 2]

        file_label.configure(text=f"Loaded File:\n{file_path}")
        messagebox.showinfo("Success", "File loaded successfully!")
    except Exception as e:
        messagebox.showerror("Error", f"Error reading file: {str(e)}")


def find_nearest_z(target):
    return unique_z_values[np.abs(np.array(unique_z_values) - target).argmin()]


def get_slice(z_target):
    return data[data["z"] == z_target]


def create_plot(column):
    if data is None:
        messagebox.showwarning(
            "Warning", "No data loaded. Please load a file first.")
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    fig.subplots_adjust(bottom=0.25, right=0.9)
    slice_data = get_slice(z_initial)

    # グリッドの作成
    x = slice_data["x"].values
    y = slice_data["y"].values
    z = slice_data[column].values

    # グリッドポイントの作成
    xi = np.linspace(x.min(), x.max(), 100)
    yi = np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)

    # 補間
    from scipy.interpolate import griddata
    zi = griddata((x, y), z, (xi, yi), method='cubic')

    # カラーマップの設定
    colormap = colormap_var.get()

    # コンタープロット
    pc = ax.pcolormesh(xi, yi, zi, cmap=colormap,
                       vmin=sali_min, vmax=sali_max, shading='auto')
    cbar = plt.colorbar(pc, label=column)

    ax.set_xlabel("X-axis(au)")
    ax.set_ylabel("Y-axis(au)")
    title = ax.set_title(
        f"J = {slice_data['jacobi constant'].values[0]:.6f}")
    jacobi_title = ax.text(
        0.5, 1.1, f"Z = {z_initial:.6f}", transform=ax.transAxes, ha='center')

    ax_slider = plt.axes([0.2, 0.1, 0.6, 0.03],
                         facecolor="lightgoldenrodyellow")
    slider = Slider(ax_slider, "Z", 0, len(unique_z_values) - 1,
                    valinit=unique_z_values.index(z_initial), valstep=1)

    def update(val):
        z_index = int(slider.val)
        z_target = unique_z_values[z_index]
        slice_data = get_slice(z_target)
        if not slice_data.empty:
            x = slice_data["x"].values
            y = slice_data["y"].values
            z = slice_data[column].values
            zi = griddata((x, y), z, (xi, yi), method='cubic')
            pc.set_array(zi.ravel())
            jacobi_title.set_text(f"Z = {z_target:.6f}")
        title.set_text(f"C_j = {slice_data['jacobi constant'].values[0]:.6f}")
        fig.canvas.draw_idle()

    slider.on_changed(update)

    coord_display = ax.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points",
                                bbox=dict(boxstyle="round", fc="w"),
                                arrowprops=dict(arrowstyle="->"))
    coord_display.set_visible(False)

    def on_hover(event):
        if event.inaxes == ax:
            x, y = event.xdata, event.ydata
            coord_display.xy = (x, y)
            coord_display.set_text(f"({x:.8f}, {y:.8f})")
            coord_display.set_visible(True)
            fig.canvas.draw_idle()
        else:
            coord_display.set_visible(False)
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", on_hover)

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

import pandas as pd
import numpy as np

# ファイルを読み込む
df = pd.read_csv("3DSALI_25_0122_2105.dat", sep="\s+", header=None, skiprows=10)

# 7列目（インデックス6）が-1でないデータのみをフィルタリング
filtered_df = df[df.iloc[:, 6] > 1]

# 7,8,9列目（インデックス6,7,8）が-1の場合、NaNに置換
filtered_df.iloc[:, 6:9] = filtered_df.iloc[:, 6:9].replace(-1, np.nan)

# 新しいファイルに書き出し
with open("pp.dat", "w") as file:

    # フィルタリングされたデータフレームを1行ずつ書き出し
    for index, row in filtered_df.iterrows():
        # 各行のデータを文字列に変換
        line = " ".join([f"{x:.12f}" if not pd.isna(x) else "nan" for x in row])
        file.write(f"{line}\n")

print("ファイルの書き出しが完了しました。")

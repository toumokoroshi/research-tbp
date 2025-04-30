set datafile separator " "  # データの区切り文字をカンマに設定

#タイトル文字の設定
set title "Trajectory sample around the Earth"
#タイトルのフォント設定
set title font"Arial,20"

set xlabel font "Arial,15"
set ylabel font "Arial,15"
set zlabel font "Arial,15"

set tics font "Arial,13"

set xlabel "X(au)" offset 0,-2 # X軸のラベル
set ylabel "Y(au)"  offset -2,0 # Y軸のラベル
set zlabel "C_j"  offset -2,0 # Z軸のラベル

set format z "%6.5f"

set size 0.9, 0.9
set origin 0.05, 0.05
set size ratio 0.5 
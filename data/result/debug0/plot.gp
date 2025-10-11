set datafile separator " "  # データの区切り文字をカンマに設定

#タイトル文字の設定
set title "Trajectory sample around the Earth"
#タイトルのフォント設定
set title font"Arial,20"

set xlabel font "Arial,15"
set ylabel font "Arial,15"
set zlabel font "Arial,15"
#ticsはメモリ文字
set tics font "Arial,15"

set xlabel "X(au)" offset 0,-2 # X軸のラベル
set ylabel "Y(au)"  offset 2,0 # Y軸のラベル
set zlabel "Z(au)"  offset 2,0 # Z軸のラベル

#set xrange [0.99:1.01]
#set yrange [-0.01:0.01]
#set zrange [-0.01:0.01]

plot 'trajectory_25_0121_1914.txt' index 0 using 2:3 w l
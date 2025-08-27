set datafile separator " "  # データの区切り文字をspaceに設定

set xlabel "X-axis(au)" offset 0,-2        # X軸のラベル
set ylabel "Y-axis(au)"  offset 5,-2       # Y軸のラベル
set zlabel "Z-axis(au)"   offset -5,0     # Z軸のラベル
#set cblabel "SALI"        # カラーバーのラベル

set xlabel font "Arial,18" 
set ylabel font "Arial,18"
set zlabel font "Arial,18"
#set cblabel font "Arial,13"  offset 8,0
set tics font "Arial,13"
#set tics font "Arial,15" tc rgb "white"
#set palette rgbformulae 33,15,10  # カラーパレットを設定（任意）
#set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#ffffff',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
#set palette defined ( 1 '#ff5555ff',1.4142 '#ff0000ff')
#set palette defined (0 '#FFFFFF',1.4142 '#ff0000ff')
set palette defined ( 0 '#FFFFFF00', 1.4142 '#FF0000FF' )

set zlabel rotate by 90
set format x "%4.3f"
set format y "%4.3f"
set format z "%4.3f"
#set colorbox        # カラーバーを表示
#set object 1 rect behind from screen 0,0 to screen 1,1 fc rgb "#333631" fillstyle solid 1.0

set view 60, 30             # 視点の角度を設定（任意）
set ticslevel 0              # Z軸の基準面を0に設定


set xrange[0.985:1.015]
set yrange[-0.015:0.015]
set zrange[-0.015:0.015]
set cbrange [2:0]

set grid
set xtics 0.99,0.005,1.01
set ytics -0.01,0.005,0.01
set ztics -0.01,0.005,0.01
set size 0.4, 0.65
set origin 0.1, 0.1
set size ratio 0.5 
unset colorbox
#splot 'output_file.txt' using 3:4:5:(column(7) == column(7) ? column(7) : 1/0) with points palette pointtype 7 pointsize 0.5
splot 'p.dat' every ::9 using 3:4:5:7 with points palette pointtype 7 pointsize 0.5
#splot '3DSALI_25_0205_1909.dat' every ::9 using 3:4:5:7 with points palette pointtype 7 pointsize 2 if (abs($7 + 1.0) > 1e-6)


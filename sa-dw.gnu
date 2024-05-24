


#set global styles
# General Palette
set palette defined ( 0 '#04233B',\
                      1 '#293E4D',\
                      2 '#20547C',\
                      3 '#77B1DE',\
                      4 '#A9C7DE',\
                      5 '#EFD5B2',\
                      6 '#EFBD79',\
                      7 '#C07F28',\
                      8 '#785F3D',\
                      9 '#5B3605')

# Reset for 3D/map plots
unset colorbox
set style line 1 pt 7 ps 0.75 lt 1 lc  rgb "#04233B" lw 2
set style line 2 pt 9 ps 0.75 lt 1 lc  rgb "#293E4D" lw 2
set style line 3 pt 5 ps 0.75 lt 1 lc  rgb "#20547C" lw 2
set style line 4 pt 11 ps 0.75 lt 1 lc rgb "#77B1DE" lw 2
set style line 5 pt 13 ps 0.75 lt 1 lc rgb "#EFBD79" lw 2
set style line 6 pt 15 ps 0.75 lt 1 lc rgb "#C07F28" lw 2
set style line 7 pt 9 ps 0.75 lt 1 lc  rgb "#785F3D" lw 2
set style line 8 pt 11 ps 0.75 lt 1 lc rgb "#5B3605" lw 2
set linetype 2 dt 2
#Set additional styles for dashed/dotted lines
set style line 100 pt 1 ps 1.2 lt 0 lc rgb "gray30" lw 2
set style line 101 pt 9 ps 1.4 lt 2 lc rgb "black" lw 2


set terminal pngcairo font "helvetica, 14"
chck_up(sl_x, sl_y, x,sl_z) = (sl_x != 1 || sl_y != 1) ? (1/0) : ((x <= 0) ? (1/0):sl_z)
chck_dw(sl_x, sl_y, x,sl_z) = (sl_x != 1 || sl_y != 1) ? (1/0) : ((x >= 0) ? (1/0):sl_z)

chck_dw_up(sl_x, sl_z, x) = (sl_z == 5 && x > 0.1) ? (sl_x) : (1/0)
chck_dw_dw(sl_x, sl_z, x) = (sl_z == 6 && x < -0.1) ? (sl_x): (1/0)

delta_S(t,m, sa) = ((t-sa*m)/sa)


SA = 1.48e7
cell_x = 0.6656

set xrange [200:300]
file= "spin-acc/0"
set output "sa-dw.png"
set ylabel "sa"
set multiplot layout 2,1

set key outside top center horizontal
set size 1,0.5
set xlabel "position (nm)"
plot file u (chck_dw_up($1,$3,$4)*cell_x):($7/SA) w p ls 3 title "sa_x",\
"" u (chck_dw_up($1,$3,$4)*cell_x):($8/SA) w p ls 4 title "sa_y",\
"" u (chck_dw_up($1,$3,$4)*cell_x):($9/SA) w p ls 5  title "sa_z"


unset title 
set size 1,0.45
set xlabel "Mn_1 SA ratio"

plot file u (chck_dw_up($1,$3,$4)*cell_x):(delta_S($7,$4, SA)) w p ls 3 title "dS_x",\
"" u (chck_dw_up($1,$3,$4)*cell_x):(delta_S($8,$5, SA)) w p ls 4 title "dS_y",\
"" u (chck_dw_up($1,$3,$4)*cell_x):(delta_S($9,$6, SA)) w p ls 5 title "dS_z",\

#"" u (chck_dw_dw($1,$3,$4)):(delta_S($7,$8,$4,$5, SA)) w p ls 6 notitle "sa_y/sa_x",\

#"" u (0.045):0 w l ls 1 dt " - " title "ab initio, dSz/dSy",\
"" u (0.29):0 w l ls 1 title 'dSx/dSy'


unset multiplot 
set auto x
set output "J_s_up.png"

set multiplot layout 2,1
#set title "spin current"
set key outside top center horizontal
set size 1,0.5
set xlabel "Mn_1 J_s (A/m^3)"
plot file u 10:(chck_up($1,$2,$4,$3)) w p ls 4 title "Js_x",\
"" u 11:(chck_up($1,$2,$4,$3)) w p ls 5 title "Js_y",\
"" u 12:(chck_up($1,$2,$4,$3)) w p ls 7 title "Js_z"

unset key
unset title 
set size 1,0.45
set xlabel "Mn_2 J_s (A/m^3)"
plot file u 10:(chck_dw($1,$2,$4,$3)) w p ls 4 notitle "sa_x",\
"" u 11:(chck_dw($1,$2,$4,$3)) w p ls 5 notitle "sa_y",\
"" u 12:(chck_dw($1,$2,$4,$3)) w p ls 7 notitle "sa_z"

unset multiplot 

set output "J_s_down.png"
set auto x
set multiplot layout 2,1
#set title "spin current"
set key outside top center horizontal
set size 1,0.5
set xlabel "Mn_1 J_s (A/m^3)"
plot file u 13:(chck_up($1,$2,$4,$3)) w p ls 4 title "Js_x",\
"" u 14:(chck_up($1,$2,$4,$3)) w p ls 5 title "Js_y",\
"" u 15:(chck_up($1,$2,$4,$3)) w p ls 7 title "Js_z"

unset key
unset title 
set size 1,0.45
set xlabel "Mn_2 J_s (A/m^3)"
plot file u 13:(chck_dw($1,$2,$4,$3)) w p ls 4 notitle "sa_x",\
"" u 14:(chck_dw($1,$2,$4,$3)) w p ls 5 notitle "sa_y",\
"" u 15:(chck_dw($1,$2,$4,$3)) w p ls 7 notitle "sa_z"

unset multiplot 
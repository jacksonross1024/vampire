


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
chck_up(sl_x, sl_y, x,sl_z) = (sl_x != 1 || sl_y != 0) ? (1/0) : ((x <= 0) ? (1/0):sl_z)
chck_dw(sl_x, sl_y, x,sl_z) = (sl_x != 1 || sl_y != 0) ? (1/0) : ((x >= 0) ? (1/0):sl_z)

chck_dw_up(sl_x, sl_z, x) = (sl_z == 5) ? (sl_x) : (1/0)
chck_dw_down(sl_x, sl_z, x) = (sl_z == 6) ? (sl_x): (1/0)

delta_S(t,m, sa) = ((t-sa*m)/sa)


SA = 1.48e7
cell_x = 1
set ytics out mirror 
set xtics out nomirror 
set ytics 0.5

fname = 4

p1 = 210
set xrange [p1-50:p1+50]

file= sprintf("spin-acc/%.0f",fname)
set output sprintf("sa-dw-sot-10e11-%.0f.png",fname)

set ylabel "sa (x 1.48e7 C/m^3)"
set multiplot layout 2,1
set title "DW STT+SOT Only"
set key outside top center horizontal
set size 1,0.5
set xlabel "position (nm)"

plot file u (chck_dw_up($1,$3,$4)*cell_x):($7/SA) w p ls 3 title "sa_x",\
"" u (chck_dw_up($1,$3,$4)*cell_x):($8/SA) w p ls 4 title "sa_y",\
"" u (chck_dw_up($1,$3,$4)*cell_x):($9/SA) w p ls 5  title "sa_z"

set ytics 0.005
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

set ytics 0.5
set xrange [p1-50:p1+50]
set output sprintf("J-sx-up-sot-10e11-%.0f.png",fname)

set ylabel "Js^x (x 10^{10} A/m^2)"
set xlabel "Position (nm)"

J_0 = 1e11
set multiplot layout 2,1
#set title "spin current"
set key outside top center horizontal
set size 1,0.5
set xlabel "Mn_1 Position (nm)"
plot file u (chck_dw_up($1,$3,$4)*cell_x):($10/J_0) w p ls 3 title "Js_x",\
"" u (chck_dw_up($1,$3,$4)*cell_x):($11/J_0) w p ls 4 title "Js_y",\
#"" u (chck_dw_up($1,$3,$4)*cell_x):($12/J_0) w p ls 5 title "Js_z",\

unset key
unset title 
set size 1,0.45
set xlabel "Mn_2 Position (nm)"
plot file u (chck_dw_down($1,$3,$4)*cell_x):($10/J_0)  w p ls 3 title "Js_x",\
"" u (chck_dw_down($1,$3,$4)*cell_x):($11/J_0)  w p ls 4 title "Js_y",\
#"" u (chck_dw_down($1,$3,$4)*cell_x):($12/J_0)  w p ls 5 title "Js_z",\

unset multiplot 

set output sprintf("J_sy_up-sot-10e11-%.0f.png",fname)

set multiplot layout 2,1
#set title "spin current"
set key outside top center horizontal
set size 1,0.5
set xlabel "Mn_1 Position (nm)"
plot file u (chck_dw_up($1,$3,$4)*cell_x):($13/J_0)  w p ls 3 title "Js_x",\
"" u (chck_dw_up($1,$3,$4)*cell_x):($14/J_0)  w p ls 4 title "Js_y",\
#"" u (chck_dw_up($1,$3,$4)*cell_x):($15/J_0)  w p ls 5 title "Js_z",\

unset key
unset title 
set size 1,0.45
set xlabel "Mn_2 Position (nm)"
plot file u (chck_dw_down($1,$3,$4)*cell_x):($13/J_0)  w p ls 3 title "Js_x",\
"" u (chck_dw_down($1,$3,$4)*cell_x):($14/J_0)  w p ls 4 title "Js_y",\
#"" u (chck_dw_down($1,$3,$4)*cell_x):($15/J_0)  w p ls 5 title "Js_z",\

unset multiplot 


set output sprintf("J_sy_down-sot-10e11-%.0f.png",fname)

set multiplot layout 2,1
#set title "spin current"
set key outside top center horizontal
set size 1,0.5
set xlabel "Mn_1 Position (nm)"
plot file u (chck_dw_up($1,$3,$4)*cell_x):($16/J_0)  w p ls 3 title "Js_x",\
"" u (chck_dw_up($1,$3,$4)*cell_x):($17/J_0)  w p ls 4 title "Js_y",\
#"" u (chck_dw_up($1,$3,$4)*cell_x):($18/J_0)  w p ls 5 title "Js_z",\

unset key
unset title 
set size 1,0.45
set xlabel "Mn_2 Position (nm)"
plot file u (chck_dw_down($1,$3,$4)*cell_x):($16/J_0)  w p ls 3 title "Js_x",\
"" u (chck_dw_down($1,$3,$4)*cell_x):($17/J_0)  w p ls 4 title "Js_y",\
#"" u (chck_dw_down($1,$3,$4)*cell_x):($18/J_0)  w p ls 5 title "Js_z",\

unset multiplot 

set output sprintf("J_sx_total-sttsot-10e11-%.0f.png",fname)

J_0 = 1e9

set ytics 10
set multiplot layout 2,1
#set title "spin current"
set key outside top center horizontal
set size 1,0.5
set xlabel "Mn_1 Position (nm)"
set title "STT + SOT"
set ylabel "Js_z (x 10^9 A/m^2)"

plot file u (chck_dw_up($1,$3,$4)*cell_x):($12/J_0)  w p ls 3 title "Js^x_z",\
"" u (chck_dw_up($1,$3,$4)*cell_x):($15/J_0)  w p ls 4 title "Js^y_z",\
"" u (chck_dw_up($1,$3,$4)*cell_x):($18/J_0)  w p ls 5 title "Js^{-y}_z",\

unset key
unset title 
set size 1,0.45
set xlabel "Mn_2 Position (nm)"
plot file u (chck_dw_down($1,$3,$4)*cell_x):($12/J_0)  w p ls 3 title "Js^x_z",\
"" u (chck_dw_down($1,$3,$4)*cell_x):($15/J_0)  w p ls 4 title "Js^y_z",\
"" u (chck_dw_down($1,$3,$4)*cell_x):($18/J_0)  w p ls 5 title "Js^{-y}_z",\

unset multiplot 


set output sprintf("Torque-stt-sot-10e11-%.0f.png",fname)

set multiplot layout 2,1
set auto y 
set ytics 3
set key outside top center horizontal
set size 1,0.5
set ylabel "Torque (mT)"
set title "STT + SOT"
set xlabel "Mn_1 Position (nm)"
plot file u (chck_dw_up($1,$3,$4)*cell_x):($27*100)  w p ls 3 title "T_x",\
"" u (chck_dw_up($1,$3,$4)*cell_x):($28*100)  w p ls 4 title "T_y",\
"" u (chck_dw_up($1,$3,$4)*cell_x):($29*100)  w p ls 5 title "T_z",\

unset key
unset title 
set size 1,0.45
set xlabel "Mn_2 Position (nm)"
plot file u (chck_dw_down($1,$3,$4)*cell_x):($27*100)  w p ls 3 title "T_x",\
"" u (chck_dw_down($1,$3,$4)*cell_x):($28*100)  w p ls 4 title "T_y",\
"" u (chck_dw_down($1,$3,$4)*cell_x):($29*100)  w p ls 5 title "T_z",\

unset multiplot 

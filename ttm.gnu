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

set style line 1 pt 7 ps 1.2 lt 1 lc  rgb "#04233B" lw 2
set style line 2 pt 9 ps 1.6 lt 1 lc  rgb "#293E4D" lw 2
set style line 3 pt 5 ps 1.2 lt 1 lc  rgb "#20547C" lw 2
set style line 4 pt 11 ps 1.6 lt 1 lc rgb "#77B1DE" lw 2
set style line 5 pt 13 ps 1.6 lt 1 lc rgb "#EFBD79" lw 2
set style line 6 pt 15 ps 1.6 lt 1 lc rgb "#C07F28" lw 2
set style line 7 pt 9 ps 1.6 lt 1 lc  rgb "#785F3D" lw 2
set style line 8 pt 11 ps 1.6 lt 1 lc rgb "#5B3605" lw 2

set linetype 2 dt 2
#Set additional styles for dashed/dotted lines
set style line 100 pt 1 ps 1.2 lt 0 lc rgb "gray30" lw 2
set style line 101 pt 9 ps 1.4 lt 2 lc rgb "black" lw 2

mb(x,t) = 0.15*sqrt(2.0/3.141592)*x*x*exp(-0.5*x*x/(t*t))/(t*t*t)
mb_e(x,t) = sqrt(2.0/3.141592)*x*x*exp(-1.0*x/t)/(t**1.5)
mb_TH(x,t) = (1.0/4.0)*sqrt(1*x/3.141592)*((1/t)**1.5)*exp(-x/t)
mb_S(x,t) = sqrt(3.141592)*(3.141592/16.0)*(x**2.0)*exp(-0.5*x/t)/(t**1.5)
FD(x) = h_1*1.0/(exp(x/TEMP_1) + 1.0)

bin(x,w) = w*floor(x/w)

kB = 1.34e-3
e  = 83.3131
n  = 97556

T0 = 400.0
t0 = T0*kB/e

s = 40.0

core = 30.86/e
transport = 80.5/e
w = 0.004

name = 0
set title "300 fs"
set ylabel "occupation"
set xlabel "e/E_f"
set boxwidth w
set key bottom left title "bin; occupation, temp"
set auto x
set auto y

#set term aqua 1
#p "Init_E_distrib" u ($0*w + core):(bin($2/e,w)) smooth freq w boxes notitle

set style fill solid 0.3 noborder
set term pngcairo font "arial, 14" size 1024,768
set output sprintf("Hist_data/%.0f.png", name)
set tics out nomirror

stats sprintf("Temp_Map0/%.0f", name) u 1:2 
N_1 = floor(STATS_records)

stats sprintf("Temp_Map1/%.0f", name) u 1:2
N_2 = floor(STATS_records)

stats sprintf("Temp_Map2/%.0f", name) u 1:2
N_3 = floor(STATS_records)

stats sprintf("Temp_Map3/%.0f", name) u 1:2
N_4 = floor(STATS_records)

stats sprintf("Temp_Map4/%.0f", name) u 1:2
N_5 = floor(STATS_records)

stats sprintf("Temp_Map5/%.0f", name) u 1:2
N_6 = floor(STATS_records)

stats sprintf("Temp_Map6/%.0f", name) u 1:2
N_7 = floor(STATS_records)

stats sprintf("Temp_Map7/%.0f", name) u 1:2
N_8 = floor(STATS_records)


h = 10

h_1 = h
TEMP_1 = 300.0
h_2 = h
TEMP_2 = 300.0
h_3 = h
TEMP_3 = 300.0
h_4 = h
TEMP_4 = 300.0
h_5 = h
TEMP_5 = 300.0
h_6 = h
TEMP_6 = 300.0
h_7 = h
TEMP_7 = 300.0
h_8 = h
TEMP_8 = 300.0

FD_1(x) = h_1/(exp((x-1.0)/TEMP_1/e) + 1.0)
FD_2(x) = h_2/(exp((x-1.0)/TEMP_2/e) + 1.0)
FD_3(x) = h_3/(exp((x-1.0)/TEMP_3/e) + 1.0)
FD_4(x) = h_4/(exp((x-1.0)/TEMP_4/e) + 1.0)
FD_5(x) = h_5/(exp((x-1.0)/TEMP_5/e) + 1.0)
FD_6(x) = h_6/(exp((x-1.0)/TEMP_6/e) + 1.0)
FD_7(x) = h_7/(exp((x-1.0)/TEMP_7/e) + 1.0)
FD_8(x) = h_8/(exp((x-1.0)/TEMP_8/e) + 1.0)

fit [core:1.03] FD_1(x) sprintf("Temp_Map0/%.0f", name) u ($1*w +core):2 via TEMP_1
fit [core:1.03] FD_2(x) sprintf("Temp_Map1/%.0f", name) u ($1*w +core):2 via TEMP_2
fit [core:1.03] FD_3(x) sprintf("Temp_Map2/%.0f", name) u ($1*w +core):2 via TEMP_3
fit [core:1.03] FD_4(x) sprintf("Temp_Map3/%.0f", name) u ($1*w +core):2 via TEMP_4
fit [core:1.03] FD_5(x) sprintf("Temp_Map4/%.0f", name) u ($1*w +core):2 via TEMP_5
fit [core:1.03] FD_6(x) sprintf("Temp_Map5/%.0f", name) u ($1*w +core):2 via TEMP_6
fit [core:1.03] FD_7(x) sprintf("Temp_Map6/%.0f", name) u ($1*w +core):2 via TEMP_7
fit [core:1.03] FD_8(x) sprintf("Temp_Map7/%.0f", name) u ($1*w +core):2 via TEMP_8

avg_T = (1.0/8.0)*(TEMP_1+TEMP_2+TEMP_3+TEMP_4+TEMP_5+TEMP_6+TEMP_7+TEMP_8)
avg_N = (1.0/8.0)*(N_1+N_2+N_3+N_4+N_5+N_6+N_7+N_8)
FD(x,t) = 1.0/(exp(x/t) + 1.0)

set yrange [0:1.5]
set xrange [0.99*core:1.1]
plot sprintf("Temp_Map0/%.0f", name) u   ($1*w + core):($2/h) smooth freq w boxes title sprintf("0; %.0f, %.0f K", N_1, TEMP_1),\
sprintf("Temp_Map1/%.0f", name) u   ($1*w + core):($2/h) smooth freq w boxes title sprintf("1; %.0f, %.0f K", N_2, TEMP_2),\
sprintf("Temp_Map2/%.0f", name) u   ($1*w + core):($2/h) smooth freq w boxes title sprintf("2; %.0f, %.0f K", N_3, TEMP_3),\
sprintf("Temp_Map3/%.0f", name) u   ($1*w + core):($2/h) smooth freq w boxes title sprintf("3; %.0f, %.0f K", N_4, TEMP_4),\
sprintf("Temp_Map4/%.0f", name) u   ($1*w + core):($2/h) smooth freq w boxes title sprintf("4; %.0f, %.0f K", N_5, TEMP_5),\
sprintf("Temp_Map5/%.0f", name) u   ($1*w + core):($2/h) smooth freq w boxes title sprintf("5; %.0f, %.0f K", N_6, TEMP_6),\
sprintf("Temp_Map6/%.0f", name) u   ($1*w + core):($2/h) smooth freq w boxes title sprintf("6; %.0f, %.0f K", N_7, TEMP_7),\
sprintf("Temp_Map7/%.0f", name) u   ($1*w + core):($2/h) smooth freq w boxes title sprintf("7; %.0f, %.0f K", N_8, TEMP_8),\
"Init_E_distrib" u ($0*w + core):(FD($0*w -1+ core, avg_T*kB/e)) w l title sprintf("avg: %.0f, %.0f K", avg_N, avg_T),\
"Init_E_distrib" u ($0*w + core):(FD($0*w -1+ core, t0)) w l lw 4 title sprintf("FEG: %.0f; %.0f K", ((30.0/4.16)**3)*4.0, T0),\
"" u (transport):($0*w) w l lw 1 notitle,\
"" u (transport + 0.004):($0*w) w l lw 1 notitle,\
"" u (0.5):($0*w) w l lw 1 notitle,\
"" u (0.524):($0*w) w l lw 1 notitle

#set print "electron_movie.txt" append
#print "file Hist_data/".name.".png"
#print "duration 0.1"

set output "TTM.png"
set auto x
set auto y

plot "mean_data.csv" u ($1*1e-3):5 w l ls 1 title "CASTLE",\
"" u ($1*1e-3):6 w l ls 2 notitle "CASTLE",\
"" u ($1*1e-3):7 w l dt 1  title "TTM",\
"" u ($1*1e-3):8 w l dt 1 notitle "TTM"



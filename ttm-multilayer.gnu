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

bin(x,w) = w*floor(x/w)

kB = 1.34e-3
e  = 108.099
n  = 9668.0

T0 = 1300.0
t0 = T0*kB/e
T1 = 300.0
t1 = T1*kB/e
s = 1e-3
core = 84.925/e
transport =  105.304/e
w = 1.0/e

name = 6000.0

set term pngcairo font "arial, 14" size 1024,768
set output "applied-field-temp.png"
set auto x
set auto y
set auto y2
set y2tics 
set tics nomirror
set key top right
set ylabel "temperature (K)"
set y2label "normalised dU"
set xlabel "time (ps)"
set title "TTM vs BTE EQ check"
stddev = 4.0*1.25*3.14*(15.0/3.53)**3

plot "mean_data.csv" u ($1*s):9 w l ls 1  title "analytic Te_1",\
"" u ($1*s):(300*sqrt($5/7.0112513386e+01)) w l ls 2 dt " . " title "global psuedo Te_1",\
"" u ($1*s):(300*sqrt($7/6.7963855344e+02)) w l ls 3 dt "-" title "local psuedo Te_1",\
"" u ($1*s):15 w l ls 4 dt " - " title "TTM Te_1",\
"" u ($1*s):11 w l ls 5  title "analytic Te_2",\
"" u ($1*s):(300*sqrt($6/1.3862321952e+02)) w l ls 6 dt " . " title "global psuedo Te_2",\
"" u ($1*s):(300*sqrt($8/1.3635270118e+03)) w l ls 7 dt "-" title "local psuedo Te_2",\
"" u ($1*s):17 w l ls 8 dt " - " title "TTM Te_2",\

#"" u ($1*s):(sqrt( ((0.46*$8+ 0.54*48.13)-(3.0*e/5.0))*4.0*e/kB/kB/3.14/3.14)) axis x1y2 w l ls 7 dt " - " title "avg KE Te",\
"" u ($1*s):(( ((0.46*$8+ 0.54*48.13)-(3.0*e/5.0))*4.0*e/kB/kB/3.14/3.14)/300.0) axis x1y2 w l ls 8 dt "  -  " title "avg KE Te constant C_v",\
#"" u ($1*s):($4/ 2.5613342888e+02) axis x1y2 w l ls 5 title "global dU",\
"" u ($1*s):($5/4.5961827195e+02) axis x1y2 w l ls 8 title "local dU",\
#"" u ($1*s):($23) axis x1y2 w p title "current density"


set ylabel "scattering count"
set y2label ""
set output "applied-field-count.png"
p "mean_data.csv" u ($1*1e-3):35 axis x1y1 w p title "electron-electron scattering < E_f",\
"" u ($1*1e-3):36 axis x1y1 w p title "electron-electron scattering > E_f",\
"" u ($1*1e-3):39 axis x1y2 w p title "electron-phonon scattering -> Te",\
"" u ($1*1e-3):40 axis x1y2 w p title "electron-phonon scattering -> Tp",\
"" u ($1*1e-3):37 axis x1y1 w p title "electron-electron scattering < E_f",\
"" u ($1*1e-3):38 axis x1y1 w p title "electron-electron scattering > E_f",\
"" u ($1*1e-3):41 axis x1y2 w p title "electron-phonon scattering -> Te",\
"" u ($1*1e-3):42 axis x1y2 w p title "electron-phonon scattering -> Tp",\


set title sprintf("%.0f fs, TTM multilayer 1", name*0.1)
set ylabel "occupation/(V * n_f /{/Symbol D}E)"
set xlabel "e/E_f"
set boxwidth w
set key top right title "Temperature"
set auto x
set auto y
set style fill solid 0.3 noborder
set output sprintf("%.0f_1.png", name)
set tics out nomirror
h = 19
FD(x,t) = 1.0/(exp(x/t) + 1.0)
set yrange [0:2]
set xrange [0.9*core:1.51]

plot sprintf("Temp_Map_1/%.0f", name) \
   u ($1*w + (0.5/e)):($4/h) smooth freq w boxes ls 1 notitle,\
"" u ($1*w + (0.5/e)):($5/h) smooth freq w boxes ls 3 notitle ,\
"" u ($1*w + (0.5/e)):($6/h) smooth freq w boxes ls 5 notitle ,\
"" u ($1*w + (0.5/e)):($7/h) smooth freq w boxes ls 7 notitle ,\
"" u ($1*w + (0.5/e)):($2/496.0) smooth freq w boxes ls 1 notitle,\
"" u ($1*w + (0.5/e)):($2/$3) w l ls 1 title "dD(e)/dt",\
"Init_E_distrib" u ($0*w + core):(FD($0*w -1+ core, t1)) w l ls 1 title sprintf("Fermi-Dirac: %.0f K", T1),\
"Init_E_distrib" u ($0*w + core):(FD($0*w -1+ core, t0)) w l ls 8 title sprintf("Fermi Dirac: %.0f K", T0),\
"" u (transport):($0*w) w l dt " - " title "global cutoff",\
#   


set title sprintf("%.0f fs, TTM multilayer 2", name*0.1)
#set output sprintf("%.0f_2.png", name)
#plot sprintf("Temp_Map_2/%.0f", name) \
   u ($1*w + (0.5/e)):($4/h) smooth freq w boxes ls 1 notitle,\
"" u ($1*w + (0.5/e)):($5/h) smooth freq w boxes ls 3 notitle ,\
"" u ($1*w + (0.5/e)):($6/h) smooth freq w boxes ls 5 notitle ,\
"" u ($1*w + (0.5/e)):($7/h) smooth freq w boxes ls 7 notitle ,\
"" u ($1*w + (0.5/e)):($2/336.0) smooth freq w boxes ls 1 notitle,\
"" u ($1*w + (0.5/e)):($2/$3) w l ls 1 title "dD(e)/dt",\
"Init_E_distrib" u ($0*w + core):(FD($0*w -1+ core, t1)) w l ls 1 title sprintf("Fermi-Dirac: %.0f K", T1),\
"Init_E_distrib" u ($0*w + core):(FD($0*w -1+ core, t0)) w l ls 8 title sprintf("Fermi Dirac: %.0f K", T0),\
"" u (transport):($0*w) w l dt " - " title "global cutoff",\

set output sprintf("E_distrib.png", name)
#plot "Init_E_vel" \
   u   (bin($2/e, 1.0/e)):(1.0/973.36) smooth freq w boxes ls 1 notitle,\
"Init_E_distrib" u ($0*w + core):(FD($0*w -1+ core, 300.0*kB/e)) w l ls 1 title sprintf("Fermi-Dirac: %.0f K", 300.0),\
"Init_E_distrib" u ($0*w + core):(FD($0*w -1+ core, t0)) w l ls 8 title sprintf("Fermi Dirac: %.0f K", T0),\
"" u (transport):($0*w) w l ls 7 notitle,\
"" u (transport + (1.0/e)):($0*w) w l ls 7 notitle,\


set auto y
set output sprintf("flux_hist_%.0f.png", name)
p sprintf("flux_hist/%.0f", name) u ($1*w + core):($2/1000.0) smooth freq w boxes title "000"

r_t(x) = A*exp(m*x)
A = 0.0035
m = 36.75
l = 1.5e3
b = 1.0
c = 1.0
set xrange [0.7:1.4]
set auto y
set auto y2
set ylabel "MFP (A)"
#set logscale y 2.718
set y2label "scattering count"
set output sprintf("relaxation_hist_%.0f.png", name)
set key top left title ""
set title "e-e scattering vs Boltzmann Transport Equilibrium"

p "Init_E_distrib" u ($1/n+1.0):(A*exp(-1*m*$1/n)) axis x1y2 w l title "{/Symbol t} = 3 fs ev^2",\
"" u ($1/n+1.0):(A*exp(-1*24.0*$1/n)) axis x1y2 w l title "{/Symbol t} = 5 fs ev^2",\
sprintf("relaxation_time_1/%.0f", name) u ($1/e):(1e-2*$2*sqrt(2.0*$1)) axis x1y1 w p title "MFP (A)",\
"" u ($1/e):($3/name) axis x1y2 w p title "scattering count",\
"Init_E_distrib" u (transport):($0*w) w l dt " - " title "global cutoff",\
sprintf("relaxation_time_2/%.0f", name) u ($1/e):(1e-2*$2*sqrt(2.0*$1)) axis x1y1 w p title "MFP (A)",\
"" u ($1/e):($3/name) axis x1y2 w p title "scattering count",\


#r_t(x) w l title sprintf("m: %.0f; A: %.0f; c: %.0f", m,A,c),\
"" u (c/2+$1/n):(3.0/(((($1/n)-c/2)**2))) w l title "1/m(x-E_f)^2",\
"" u (c/2+$1/n):(40+A-A*exp(-1*l*(($1/n)-c/2)**2)) axis x1y1 w l title "Theoretical Approximation K_{ee} = 3 fs eV^2",\
"" u (c/2+$1/n):(40+A-A*exp(-1*2e3*(($1/n)-c/2)**2)) axis x1y1 w l title "Simulation Constant K_{ee} = 5 fs eV^2",\


set output "delta_occupation.png"
set xrange [0.95:1.05]
set auto y
set title "change in occupation"
set ylabel "occupation"
set xlabel "e/E_f"
set key title "time step"
set key top right

plot for [i=40000:70000:5000] 'Temp_Map_1/'.i.'' u ($1*w + 0.75/e):($2) w l ls (i-40000.0)/5000.0 title sprintf("%f",i),\
"Init_E_distrib" u ($1*w + core):(47*(FD($0*w -1+ core, t1))) w l ls 8 dt " - " title sprintf("Fermi Dirac: %.0f K", T1),\


set output "ep-analytic.png"
set ylabel "analytic"
set xlabel "discrete"
set xrange [-0.5:1]
set title "discrete vs analytic ep scattering"
set key notitle top left
set y2label ""
#p "ee_scattering" u 2:8 w p title "forward rate",\
"" u 3:9 w p title "reverse rate",\
"" u ($2+$3):($8-$9) w p title "total rate",\
"Init_E_distrib" u ($0/n-0.25):($0/n-0.25) w l notitle,\
"ee_scattering" u 5:10 w p title "e-occupation",\
"" u 6:11 w p title "f-e-occupation",\
"" u 7:12 w p title "r-e-occupation"


set output "ee-analytic.png"
set ylabel "analytic"
set xlabel "discrete"
set auto x
set auto x2
set auto y2
set title "discrete vs analytic ee scattering"
set key notitle top left
set y2label ""
set x2label ""
set x2tics nomirror

#p "ee_scattering" \
u ($5):4 w p title "H(k)_{sr}",\
"" u ($6):4 axis x2y2 w p title "H(e)_{sr}",\


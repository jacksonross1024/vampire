

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


set term pngcairo font "helvetica, 14" size 600, 600
set xlabel "Time (ps)"
set ylabel "Depth (nm)"
set ytics out nomirror
set xtics out nomirror
set xrange [-0.2:4]
cell = 0.608

tp = 0.085 # ps
a = 25.0 #penetration, nm
h = 60.8 #nm
time = 1000 # ps
set yrange [0:h]
file = "vertical_temperature_profile.dat"
set style fill solid border
gaus(t,z) = exp(-z/a)* exp(-4*log(2)*((t-2.0*tp)/tp)**2.0)
set colorbox 
set cbrange [0:]
tsi = 1000 #tsi
s = 1e-3*tsi #convert tsi to ps
o = 0.0  #offset

set cbrange [5:13]
set xrange [0:time]
set cblabel "T_p (K)"

set logscale x 10 
higher_temp(Te, Tp) = Tp

set output "TTMp-0.1Sink-debeye-model-2nd-der-reduced-lr.png"
set title "12nm h-BN; 4L Moire CrI__3; 12nm h-BN; 30 nm SiO2"
plot for [i=0:h/cell:1] file u ($1*s -o):(i*cell):(s/2.0):(cell/2.0):(higher_temp( column(2+(2*i)), column(3+(2*i)))) w boxxy lc palette notitle,\
(14.007) w l dt " - " lc "white" notitle,\
(16.443) w l dt " - " lc "white" notitle,\
(30.45) w l dt " - " lc "white" notitle


set output "TTMe-0.1Sink-debeye-model-2nd-der-reduced-lr.png"

higher_temp(Te, Tp) = Te 
set title "12nm h-BN; 4L Moire CrI__3; 12nm h-BN; 30 nm SiO2"
set cblabel "T_e (K)"
plot for [i=0:h/cell:1] file u ($1*s -o):(i*cell):(s/2.0):(cell/2.0):(higher_temp( column(2+(2*i)), column(3+(2*i)))) w boxxy lc palette notitle,\
(14.007) w l dt " - " lc "white" notitle,\
(16.443) w l dt " - " lc "white" notitle,\
(30.45) w l dt " - " lc "white" notitle

set yrange [10:20]
higher_temp(Te, Tp) = Tp
set output "TTMp-0.1Sink-debeye-model-2nd-der-reduced-interface-hr-lr.png"
set title "12nm h-BN; 4L Moire CrI__3; 12nm h-BN; 30 nm SiO2"
plot for [i=0:h/cell:1] file u ($1*s -o):(i*cell):(s/2.0):(cell/2.0):(higher_temp( column(2+(2*i)), column(3+(2*i)))) w boxxy lc palette notitle,\
(14.007) w l dt " - " lc "white" notitle,\
(16.443) w l dt " - " lc "white" notitle,\
(30.45) w l dt " - " lc "white" notitle

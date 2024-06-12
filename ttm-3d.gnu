
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


set term pngcairo font "helvetica, 14"
set xlabel "Time (ps)"
set ylabel "Depth (nm)"
set ytics out nomirror
set xtics out nomirror
set xrange [-1:0.5]
cell = 0.142316667 
set yrange [0:50]
tp = 0.05 
a = 30.0

gaus(t,z) = exp(-z/a)* exp(-4*log(2)*((t-2.0*tp)/tp)**2.0)
set colorbox 

#set output "TTM-gaussian-attentuation-lot.png"
set clabel "Field (mT)"
#plot for [i=0:352:10] "vertical_temperature_profile.dat" u ($1*1e-3-2):(50-i*cell):(3.67*gaus($1*1e-3-2, 50-i*cell)) w boxes lc palette notitle

set xrange [0:5]
set clabel "Temperature (K)"
set output "TTMe-gaussian-attenuation.png"
plot for [i=0:352:10] "vertical_temperature_profile.dat"  u ($1*1e-3 -2):(50-i*cell):2+(2*i) w boxes lc palette notitle


set output "TTMp-gaussian-attenuation.png"
set xrange [0:25]
plot for [i=0:352:10] "vertical_temperature_profile.dat" u ($1*1e-3 -2):(50-i*cell):3+(2*i) w boxes lc palette notitle

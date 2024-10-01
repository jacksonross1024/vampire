


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
set xrange [-0.2:4]
cell = 0.142316667 

tp = 1.0 
a = 30.0
h = 10.6
time = 12
set yrange [0:h]
file = "vertical_temperature_profile.dat"
set style fill solid 
gaus(t,z) = exp(-z/a)* exp(-4*log(2)*((t-1.5*tp)/tp)**2.0)
set colorbox 
set cbrange [0:]
s = 1e-2
o = 0.2 
set xrange [-0.2:12]
set output "TTM-gaussian-attentuation-lot.png"
set cblabel "Intensity (1/I_0)"
plot for [i=0:h/cell:1] file u ($1*s-o):(h-i*cell):(s):(cell):(gaus($1*s-o, h-i*cell)) w boxxy lc palette notitle

set cbrange [300:]
set xrange [0:time]
set cblabel "Temperature (K)"
set output "TTMe-gaussian-attenuation.png"
plot for [i=0:h/cell:1] file  u ($1*s -o):(h-i*cell):(s):(cell):2+(2*i) w boxxy lc palette notitle


set output "TTMp-gaussian-attenuation.png"
set xrange [0:time]
plot for [i=0:h/cell:1] file u ($1*s -o):(h-i*cell):(s):(cell):3+(2*i) w boxxy lc palette notitle


set output "Mag-gaussian-attenuation.png"
set cblabel "Magnetisation (m_x/m_s)"
set xrange [0:time]
set cbrange [-1:1]
m(i) = (5 + 4*i)*(2 + 4*i)

file_1 = "microcells.txt"
plot for [i=0:h/cell:2] file_1  u ($1*s -o):(h-2*i*cell):(s):(cell*2):6+i*4 w boxxy lc palette notitle

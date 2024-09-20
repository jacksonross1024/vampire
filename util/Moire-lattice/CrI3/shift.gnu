
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

set palette defined (  0 'blue', 1 'white',  2 'red')

set term pngcairo font "helvetica, 14"

set xlabel "x (A)"
set ylabel "y (A)"
set ytics out nomirror
set xtics out nomirror



a0x = 6.93
a1x = -3.465
a1y = 6.002

set style fill solid noborder


set xrange [-1:1]
set yrange [-1:1]


set colorbox 
set auto cb
unset colorbox 
# 6.930   0.000   0.000
#-3.465   6.002   0.000
# 0.000   0.000  25.000

set output "Jinter.png"


set xrange [-a0x:a0x]
set yrange [-a0x:a0x]

set multiplot 
set origin 0,0.4
set size 0.5,0.6
set key title "DFT"
bounds(x,y,z) = (sqrt(x*x + y*y) > 6.93) ? (0.0) : (z)
plot "Jintra1_interpolation_out.txt" u (a0x*($1*0.01-1)):(a0x*($2*0.01-1)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$3)) w boxxy palette notitle

set origin 0.0,-0.1
set size 0.5,0.6
set key title "Interpolation"
plot "files/interpolated_array" u (a0x*($1*0.02-1)):(a0x*($2*0.02-1)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.02-1),a0x*($2*0.02-1),$3)) w boxxy palette notitle
unset multiplot 


set output "Dinter.png"
set multiplot 
set origin 0,0.5
set size 0.5,0.5
set xrange [-a0x:a0x]
set yrange [-a0x:a0x]


set key title "DFT"
plot "Maps/Interpolated_Dij_Inter" u (a0x*$1*0.01-a0x):(a1y*($2*0.01)-a1y):(a0x*0.01):(a1y*0.01):($5) w boxxy palette notitle

set origin 0,0.0
set size 0.5,0.5
set key title "Interpolation"
plot "files/interpolated_Dij_inter" u (a0x*$1*0.01-a0x):(a0x*$2*0.01-a0x):(a0x*0.01):(a1y*0.01):($5) w boxxy palette notitle
unset multiplot 

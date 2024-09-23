
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

set term pngcairo font "helvetica, 14" size 900,900

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

unset colorbox 
# 6.930   0.000   0.000
#-3.465   6.002   0.000
# 0.000   0.000  25.000

set output "Dintera-rotation.png"

set xrange [-a0x:a0x]
set yrange [-a0x:a0x]
set angles degrees
set multiplot layout 3,3
rot = 90.0
rot2 = rot + 180.0
set key title "DFT"
bounds(x,y,z) = (sqrt(x*x + y*y) > 6.93) ? (0.0) : (z)
set cbrange [-0.1:0.1]
plot "files/criteria.txt" u (a0x*$1):(a1y*($2)):(a0x*0.05):(a1y*0.0433):($9) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1):(a1y*($2)):(a0x*0.05):(a1y*0.0433):($10) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1):(a1y*($2)):(a0x*0.05):(a1y*0.0433):($11) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1*cos(rot)-a0x*$2*sin(rot)):(a1y*($2)*cos(rot)+a1y*sin(rot)*$1):(a0x*0.05):(a1y*0.0433):($9*cos(rot)-$10*sin(rot)) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1*cos(rot)-a0x*$2*sin(rot)):(a1y*($2)*cos(rot)+a1y*sin(rot)*$1):(a0x*0.05):(a1y*0.0433):($10*cos(rot)+$9*sin(rot)) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1*cos(rot)-a0x*$2*sin(rot)):(a1y*($2)*cos(rot)+a1y*sin(rot)*$1):(a0x*0.05):(a1y*0.0433):($11) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1*cos(rot2)-a0x*$2*sin(rot2)):(a1y*($2)*cos(rot2)+a1y*sin(rot2)*$1):(a0x*0.05):(a1y*0.0433):($9*cos(rot2)-$10*sin(rot2)) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1*cos(rot2)-a0x*$2*sin(rot2)):(a1y*($2)*cos(rot2)+a1y*sin(rot2)*$1):(a0x*0.05):(a1y*0.0433):($10*cos(rot2)+$9*sin(rot2)) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1*cos(rot2)-a0x*$2*sin(rot2)):(a1y*($2)*cos(rot2)+a1y*sin(rot2)*$1):(a0x*0.05):(a1y*0.0433):($11) w boxxy palette notitle


#set key title "Interpolation"
#plot "files/interpolated_array" u (a0x*($1*0.02-1)):(a1y*($2*0.02-1)):(a0x*0.02):(a1y*0.02):(bounds(a0x*($1*0.02-1),a1y*($2*0.02-1),$3)) w boxxy palette notitle
unset multiplot 


set output "Dinter.png"
set multiplot layout 3,3

set xrange [-a0x:a0x]
set yrange [-a0x:a0x]


set key title "DFT"
plot "files/interpolated_Dij_inter" u (a0x*$1*0.01-a0x):(a0x*($2*0.01)-a0x):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$3)) w boxxy palette notitle

plot "files/interpolated_Dij_inter" u (a0x*$1*0.01-a0x):(a0x*($2*0.01)-a0x):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$4)) w boxxy palette notitle

plot "files/interpolated_Dij_inter" u (a0x*$1*0.01-a0x):(a0x*($2*0.01)-a0x):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$5)) w boxxy palette notitle

set key title "Interpolation"
plot "Maps/Interpolated-Data/Interpolated_Dij_Inter.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$3)) w boxxy palette notitle

plot "Maps/Interpolated-Data/Interpolated_Dij_Inter.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$4)) w boxxy palette notitle

plot "Maps/Interpolated-Data/Interpolated_Dij_Inter.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$5)) w boxxy palette notitle

plot "files/Dinter.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle

plot "files/Dinter.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$4)) w boxxy palette notitle

plot "files/Dinter.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$5)) w boxxy palette notitle


unset multiplot 

set output "Dintra.png"
set multiplot layout 3,3

set xrange [-a0x:a0x]
set yrange [-a0x:a0x]


set key title "DFT"
plot "files/interpolated_Dij_intra" u (a0x*$1*0.01-a0x):(a0x*($2*0.01)-a0x):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$3)) w boxxy palette notitle

plot "files/interpolated_Dij_intra" u (a0x*$1*0.01-a0x):(a0x*($2*0.01)-a0x):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$4)) w boxxy palette notitle

plot "files/interpolated_Dij_intra" u (a0x*$1*0.01-a0x):(a0x*($2*0.01)-a0x):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$5)) w boxxy palette notitle

set key title "Interpolation"
plot "Maps/Interpolated-Data/Interpolated_1st_Dij_Intra.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$3)) w boxxy palette notitle

plot "Maps/Interpolated-Data/Interpolated_1st_Dij_Intra.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$4)) w boxxy palette notitle

plot "Maps/Interpolated-Data/Interpolated_1st_Dij_Intra.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$5)) w boxxy palette notitle

plot "files/Dintra-correct.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle

plot "files/Dintra-correct.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$4)) w boxxy palette notitle

plot "files/Dintra-correct.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$5)) w boxxy palette notitle


unset multiplot 


set output "Dintra2.png"
set multiplot layout 3,3

set xrange [-a0x:a0x]
set yrange [-a0x:a0x]

set cbrange [-0.20:0.20]
set key title "DFT"
plot "files/interpolated_Dij_intra" u (a0x*$1*0.01-a0x):(a0x*($2*0.01)-a0x):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$3)) w boxxy palette notitle

plot "files/interpolated_Dij_intra" u (a0x*$1*0.01-a0x):(a0x*($2*0.01)-a0x):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$4)) w boxxy palette notitle

plot "files/interpolated_Dij_intra" u (a0x*$1*0.01-a0x):(a0x*($2*0.01)-a0x):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$5)) w boxxy palette notitle

set key title "Interpolation"
plot "Maps/Interpolated-Data/Interpolated_2nd_Dij_Intra.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$3)) w boxxy palette notitle

plot "Maps/Interpolated-Data/Interpolated_2nd_Dij_Intra.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$4)) w boxxy palette notitle

plot "Maps/Interpolated-Data/Interpolated_2nd_Dij_Intra.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$5)) w boxxy palette notitle

plot "files/Dintra-correct.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle

plot "files/Dintra-correct.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$4)) w boxxy palette notitle

plot "files/Dintra-correct.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$5)) w boxxy palette notitle


unset multiplot 


set cbrange [1.0:3.0]
set output "Jintra.png"
set multiplot layout 3,3

set xrange [-a0x:a0x]
set yrange [-a0x:a0x]


set key title "DFT"
plot "files/interpolated_J1_intra_AA.txt" u (a0x*$1*0.02-a0x):(a0x*($2*0.02)-a0x):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.02-1),a0x*($2*0.02-1),$3)) w boxxy palette notitle

plot "files/interpolated_J1_intra_AA.txt" u (a0x*$1*0.02-a0x):(a0x*($2*0.02)-a0x):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.02-1),a0x*($2*0.02-1),$3)) w boxxy palette notitle

plot "files/interpolated_J1_intra_AA.txt" u (a0x*$1*0.02-a0x):(a0x*($2*0.02)-a0x):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.02-1),a0x*($2*0.02-1),$3)) w boxxy palette notitle

set key title "Interpolation"
plot "Maps/Interpolated-Data/Interpolated_J1_Intra.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$3)) w boxxy palette notitle

plot "Maps/Interpolated-Data/Interpolated_J1_Intra.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$3)) w boxxy palette notitle

plot "Maps/Interpolated-Data/Interpolated_J1_Intra.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$3)) w boxxy palette notitle

plot "files/Jintra_1N.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle

plot "files/Jintra_1N.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle

plot "files/Jintra_1N.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle


unset multiplot 

set cbrange [0.30:0.90]
set output "Jintra2.png"
set multiplot layout 3,3

set xrange [-a0x:a0x]
set yrange [-a0x:a0x]


set key title "DFT"
plot "files/interpolated_J2_intra_AA.txt" u (a0x*$1*0.02-a0x):(a0x*($2*0.02)-a0x):(a0x*0.02):(a0x*0.02):(bounds(a0x*($1*0.02-1),a0x*($2*0.02-1),$3)) w boxxy palette notitle

plot "files/interpolated_J2_intra_AA.txt" u (a0x*$1*0.02-a0x):(a0x*($2*0.02)-a0x):(a0x*0.02):(a0x*0.02):(bounds(a0x*($1*0.02-1),a0x*($2*0.02-1),$3)) w boxxy palette notitle

plot "files/interpolated_J2_intra_AA.txt" u (a0x*$1*0.02-a0x):(a0x*($2*0.02)-a0x):(a0x*0.02):(a0x*0.02):(bounds(a0x*($1*0.02-1),a0x*($2*0.02-1),$3)) w boxxy palette notitle

set key title "Interpolation"
plot "Maps/Interpolated-Data/Interpolated_J2_Intra.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$3)) w boxxy palette notitle

plot "Maps/Interpolated-Data/Interpolated_J2_Intra.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$3)) w boxxy palette notitle

plot "Maps/Interpolated-Data/Interpolated_J2_Intra.txt" u (a0x*(0.01*$2)):(a0x*(0.01*$1)):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$3)) w boxxy palette notitle

plot "files/Jintra_2N.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle

plot "files/Jintra_2N.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle

plot "files/Jintra_2N.csv" u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle


unset multiplot 

set auto y
set auto x 

set ylabel "nm"
set xlabel "nm"
set output "mag.png"
set auto cb
unset key 

plot "cells-00000000.txt" u ($1*0.1):($2*0.1):(1.0):(1.0):($22*$24) w boxxy palette notitle 
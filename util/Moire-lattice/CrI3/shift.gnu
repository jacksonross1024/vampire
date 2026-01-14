
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
set output "MvsT.png"
set ylabel "M/M_s"
set xlabel "Temperature (K)"

set xrange [0:60]
set yrange [-0.1:1.1]
set multiplot 

beta = 0.5
T_C = 35
set ytics out nomirror 
set xtics out nomirror 

mag(t, T_C, beta) = (t>T_C) ? (0) : (1-t/T_C)**beta

#fit [0:T_C] mag(x, T_C, beta) "output-0.0-3J-DMI-mvsT" u 1:9 via T_C, beta

#fit [100:T_c_1] mvsT(x,a_1,b_1,T_c_1) "output-18nm3" u 1:(($4+$5)*0.5)via T_c_1

set key  top right reverse title "0.0^0"
plot "output-0.0-1J-DMI-mvsT" u 1:9 w p ls 2 title "1NN",\
"output-0.0-2J-DMI-mvsT" u 1:9 w p ls 3 title "2NN",\
"output-0.0-3J-DMI-mvsT" u 1:9 w p ls 4 title "3NN",\
mag(x, T_C, beta) w l ls 5 title sprintf("T_C: %.0f; beta: %.3f", T_C, beta)

unset multiplot

set palette defined (  0 'blue', 1 'white',  2 'red')


set term pngcairo font "helvetica, 14" size 900,900




unset colorbox 
# 6.930   0.000   0.000
#-3.465   6.002   0.000
# 0.000   0.000  25.000

set output "Dintera-rotation.png"

set xrange [-a0x:a0x]
set yrange [-a0x:a0x]
set angles degrees
set multiplot layout 3,3
rot = 120
rot2 = 90.0
set key title "DFT"
bounds(x,y,z) = (sqrt(x*x + y*y) > 6.93) ? (0.0) : (z)
set cbrange [-0.1:0.1]
plot "files/criteria.txt" u (a0x*$1):(a1y*($2)):(a0x*0.05):(a1y*0.0433):($12*cos(rot2)-$13*sin(rot2)) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1):(a1y*($2)):(a0x*0.05):(a1y*0.0433):($13*cos(rot2)+$12*sin(rot2)) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1):(a1y*($2)):(a0x*0.05):(a1y*0.0433):($14) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1*cos(rot)-a0x*$2*sin(rot)):(a1y*($2)*cos(rot)+a1y*sin(rot)*$1):(a0x*0.05):(a1y*0.0433):(cos(rot2)*($9*cos(rot)-$10*sin(rot))-sin(rot2)*($10*cos(rot)+$9*sin(rot))) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1*cos(rot)-a0x*$2*sin(rot)):(a1y*($2)*cos(rot)+a1y*sin(rot)*$1):(a0x*0.05):(a1y*0.0433):(cos(rot2)*($10*cos(rot)+$9*sin(rot))+sin(rot2)*($9*cos(rot)-$10*sin(rot))) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1*cos(rot)-a0x*$2*sin(rot)):(a1y*($2)*cos(rot)+a1y*sin(rot)*$1):(a0x*0.05):(a1y*0.0433):($11) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1*cos(rot2)-a0x*$2*sin(rot2)):(a1y*($2)*cos(rot2)+a1y*sin(rot2)*$1):(a0x*0.05):(a1y*0.0433):($12) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1*cos(rot2)-a0x*$2*sin(rot2)):(a1y*($2)*cos(rot2)+a1y*sin(rot2)*$1):(a0x*0.05):(a1y*0.0433):($13) w boxxy palette notitle

plot "files/criteria.txt" u (a0x*$1*cos(rot2)-a0x*$2*sin(rot2)):(a1y*($2)*cos(rot2)+a1y*sin(rot2)*$1):(a0x*0.05):(a1y*0.0433):($14) w boxxy palette notitle


#set key title "Interpolation"
#plot "files/interpolated_array" u (a0x*($1*0.02-1)):(a1y*($2*0.02-1)):(a0x*0.02):(a1y*0.02):(bounds(a0x*($1*0.02-1),a1y*($2*0.02-1),$3)) w boxxy palette notitle
unset multiplot 


set output "Dinter.png"
set multiplot layout 3,3

set xrange [-1.1:1.1]
set yrange [-1.1:1.1]


unset key 
file1 = "files/Interpolated_Dij_Inter"

plot file1 u ($2*0.01):(($1*0.01)*0.866):(0.01):(0.01*0.866):(bounds(($1*0.01),($2*0.01*0.866),$3)) w boxxy palette notitle

plot file1 u ($2*0.01):(($1*0.01)*0.866):(0.01):(0.01*0.866):(bounds(($1*0.01),($2*0.01*0.866),$4)) w boxxy palette notitle

plot file1 u ($2*0.01):(($1*0.01)*0.866):(0.01):(0.01*0.866):(bounds(($1*0.01),($2*0.01*0.866),$5)) w boxxy palette notitle


file2 = "files/Interpolated_Dij_Inter_out.txt"
plot file2  u ((0.01*$1)-1):(0.866*(0.01*$2)-0.866):(0.01):(0.866*0.01):(bounds(($1*0.01-1),0.866*($2*0.01-1),$3)) w boxxy palette notitle

plot file2  u ((0.01*$1)-1):(0.866*(0.01*$2)-0.866):(0.01):(0.866*0.01):(bounds(($1*0.01-1),0.866*($2*0.01-1),$4)) w boxxy palette notitle

plot file2  u ((0.01*$1)-1):(0.866*(0.01*$2)-0.866):(0.01):(0.866*0.01):(bounds(($1*0.01-1),0.866*($2*0.01-1),$5)) w boxxy palette notitle


file3 = "files/Dinter.csv"
plot file3 u ($1):($2):(0.05):(0.0866):(bounds(($1),($2),$3)) w boxxy palette notitle

plot file3 u ($1):($2):(0.05):(0.0866):(bounds(($1),($2),$4)) w boxxy palette notitle 

plot file3 u ($1):($2):(0.05):(0.0866):(bounds(($1),($2),$5)) w boxxy palette notitle


unset multiplot 


set output "Dintra.png"
set multiplot layout 3,3

set xrange [-a0x:a0x]
set yrange [-a0x:a0x]

unset key 
file1 = "files/Interpolated_1st_Dij_Intra"

plot file1 u (a0x*$2*0.01):(a1y*($1*0.01)):(a0x*0.01):(a1y*0.01):(bounds(($1*0.01),a1y*($2*0.01),$3)) w boxxy palette notitle

plot file1 u (a0x*$2*0.01):(a1y*($1*0.01)):(a0x*0.01):(a1y*0.01):(bounds(a0x*($1*0.01),a1y*($2*0.01),$4)) w boxxy palette notitle

plot file1 u (a0x*$2*0.01):(a1y*($1*0.01)):(a0x*0.01):(a1y*0.01):(bounds(a0x*($1*0.01),a1y*($2*0.01),$5)) w boxxy palette notitle


file2 = "files/Interpolated_1st_Dij_Intra_out.txt"
plot file2  u (a0x*(0.01*$1)-a0x):(a1y*(0.01*$2)-a1y):(0.01*a0x):(a1y*0.01):($3) w boxxy palette notitle

plot file2 u (a0x*(0.01*$1)-a0x):(a1y*(0.01*$2)-a1y):(0.01*a0x):(a1y*0.01):(bounds(a0x*($1*0.01-1),a1y*($2*0.01-1),$4)) w boxxy palette notitle

plot file2 u (a0x*(0.01*$1)-a0x):(a1y*(0.01*$2)-a1y):(0.01*a0x):(a1y*0.01):(bounds(a0x*($1*0.01-1),a1y*($2*0.01-1),$5)) w boxxy palette notitle


file3 = "files/Dintra-correct.csv"
plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle

plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$4)) w boxxy palette notitle

plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$5)) w boxxy palette notitle

unset multiplot 


set output "Dintra2.png"
set multiplot layout 3,3

set xrange [-a0x:a0x]
set yrange [-a0x:a0x]

set cbrange [-0.20:0.20]

unset key 
file1 = "files/Interpolated_2nd_Dij_Intra"

plot file1 u (a0x*$1*0.01):(a0x*($2*0.01)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$3)) w boxxy palette notitle

plot file1 u (a0x*$1*0.01):(a0x*($2*0.01)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$4)) w boxxy palette notitle

plot file1 u (a0x*$1*0.01):(a0x*($2*0.01)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$5)) w boxxy palette notitle


file2 = "files/Interpolated_2nd_Dij_Intra_out.txt"
plot file2  u (a0x*(0.01*$1)-a0x):(a0x*(0.01*$2)-a0x):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$3)) w boxxy palette notitle

plot file2 u (a0x*(0.01*$1)-a0x):(a0x*(0.01*$2)-a0x):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$4)) w boxxy palette notitle

plot file2 u (a0x*(0.01*$1)-a0x):(a0x*(0.01*$2)-a0x):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$5)) w boxxy palette notitle


file3 = "files/Dintra.csv"
plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle

plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$4)) w boxxy palette notitle

plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$5)) w boxxy palette notitle

unset multiplot 


set cbrange [1.0:3.0]
set output "Jintra.png"
set multiplot layout 3,3

set xrange [-a0x:a0x]
set yrange [-a0x:a0x]


unset key 
file1 = "files/Interpolated_J1_Intra"

plot file1 u (a0x*$1*0.01):(a0x*($2*0.01)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$3)) w boxxy palette notitle

plot file1 u (a0x*$1*0.01):(a0x*($2*0.01)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$4)) w boxxy palette notitle

plot file1 u (a0x*$1*0.01):(a0x*($2*0.01)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$5)) w boxxy palette notitle


file2 = "files/Interpolated_J1_Intra_out.txt"
plot file2  u (a0x*(0.01*$1)-a0x):(a0x*(0.01*$2)-a0x):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$3)) w boxxy palette notitle

plot file2 u (a0x*(0.01*$1)-a0x):(a0x*(0.01*$2)-a0x):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$4)) w boxxy palette notitle

plot file2 u (a0x*(0.01*$1)-a0x):(a0x*(0.01*$2)-a0x):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$5)) w boxxy palette notitle


file3 = "files/Jintra_1N.csv"
plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle

plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$4)) w boxxy palette notitle

plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$5)) w boxxy palette notitle

unset multiplot 


set cbrange [0.30:0.90]
set output "Jintra2.png"
set multiplot layout 3,3

set xrange [-a0x:a0x]
set yrange [-a0x:a0x]


unset key 
file1 = "files/Interpolated_J2_Intra"

plot file1 u (a0x*$1*0.01):(a0x*($2*0.01)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$3)) w boxxy palette notitle

plot file1 u (a0x*$1*0.01):(a0x*($2*0.01)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$4)) w boxxy palette notitle

plot file1 u (a0x*$1*0.01):(a0x*($2*0.01)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$5)) w boxxy palette notitle


file2 = "files/Interpolated_J2_Intra_out.txt"
plot file2  u (a0x*(0.01*$1)-a0x):(a0x*(0.01*$2)-a0x):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$3)) w boxxy palette notitle

plot file2 u (a0x*(0.01*$1)-a0x):(a0x*(0.01*$2)-a0x):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$4)) w boxxy palette notitle

plot file2 u (a0x*(0.01*$1)-a0x):(a0x*(0.01*$2)-a0x):(0.01*a0x):(a0x*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$5)) w boxxy palette notitle


file3 = "files/Jintra_2N.csv"
plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$3)) w boxxy palette notitle

plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$4)) w boxxy palette notitle

plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$5)) w boxxy palette notitle


unset multiplot 

set cbrange [-0.5:0.5]
set output "Jinter.png"
set multiplot layout 3,3


set xrange [-1.1:1.1]
set yrange [-1.1:1.1]
unset key 
file1 = "files/Interpolated_J_Inter"

plot file1 u ($2*0.01):(($1*0.01)*0.866):(0.01):(0.01*0.866):(bounds(($1*0.01),($2*0.01),$3)) w boxxy palette notitle

plot file1 u (a0x*$1*0.01):(a0x*($2*0.01)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$4)) w boxxy palette notitle

plot file1 u (a0x*$1*0.01):(a0x*($2*0.01)):(a0x*0.01):(a0x*0.01):(bounds(a0x*($1*0.01),a0x*($2*0.01),$5)) w boxxy palette notitle


file2 = "files/Interpolated_J_Inter_out.txt"
plot file2  u ($1*0.01-1):(($2*0.01-1)*0.866):(0.01):(0.01*0.866):(bounds(($1*0.01-1),($2*0.01-1),$3))  w boxxy palette notitle

plot file2 u (a0x*(0.01*$1)-a0x):(a0x*(0.01*$2)-a0x):(0.01*a0x):(a1y*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$4)) w boxxy palette notitle

plot file2 u (a0x*(0.01*$1)-a0x):(a0x*(0.01*$2)-a0x):(0.01*a0x):(a1y*0.01):(bounds(a0x*($1*0.01-1),a0x*($2*0.01-1),$5)) w boxxy palette notitle


file3 = "files/Jinter.csv"
plot file3 u ($1):($2):(0.05):(0.0866):(bounds(($1),($2),$3)) w boxxy palette notitle

plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$4)) w boxxy palette notitle

plot file3 u (a0x*$1):(a0x*$2):(0.05*a0x):(a0x*0.0866):(bounds(a0x*($1),a0x*($2),$5)) w boxxy palette notitle


unset multiplot 

set auto y
set auto x 

set ylabel "nm"
set xlabel "nm"
set output "mag.png"
set auto cb
unset key 

plot "cells-00000000.txt" u ($1*0.1):($2*0.1):(1.0):(1.0):($22*$24) w boxxy palette notitle 
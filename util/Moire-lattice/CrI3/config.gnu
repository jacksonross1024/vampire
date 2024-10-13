
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


set term pngcairo size 600,600



set ylabel "y position (nm)"
set xlabel "x position (nm)"

set ytics 50 out 
set xtics 50 out 
set mytics 5 
set mxtics 5

set colorbox

unset colorbox
set palette defined (  0 'blue', 1 'white',  2 'red')


set colorbox

set xrange [2:198]
set yrange [2:198]

file = "config_energy-0.5-2NN-DMI"
set output sprintf("%s.png", file)
set multiplot layout 2,3

set colorbox
set title "bottom layer, D_x"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

unset colorbox 
unset ylabel 
unset xlabel 

set title "bottom layer, D_y"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$3) w boxxy palette notitle 

set title "bottom layer, D_z"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$3) w boxxy palette notitle 

set ylabel "y position (nm)"
set xlabel "x position (nm)"
set colorbox
set title "top layer, D_x"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($8/$4) w boxxy palette notitle  

unset colorbox 
unset ylabel 
unset xlabel 
set title "top layer, D_y"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$4) w boxxy palette notitle  

set title "top layer, D_z"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$4) w boxxy palette notitle  

unset multiplot 

file = "config_energy-1.1-2NN-DMI"
set output sprintf("%s.png", file)
set multiplot layout 2,3

set colorbox
set title "bottom layer, D_x"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

unset colorbox 
unset ylabel 
unset xlabel 

set title "bottom layer, D_y"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$3) w boxxy palette notitle 

set title "bottom layer, D_z"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$3) w boxxy palette notitle 

set ylabel "y position (nm)"
set xlabel "x position (nm)"
set colorbox
set title "top layer, D_x"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($8/$4) w boxxy palette notitle  

unset colorbox 
unset ylabel 
unset xlabel 
set title "top layer, D_y"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$4) w boxxy palette notitle  

set title "top layer, D_z"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$4) w boxxy palette notitle  

unset multiplot 

file = "config_energy-0.5-3NN-DMI"
set output sprintf("%s.png", file)
set multiplot layout 2,3

set colorbox
set title "bottom layer, D_x"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

unset colorbox 
unset ylabel 
unset xlabel 

set title "bottom layer, D_y"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$3) w boxxy palette notitle 

set title "bottom layer, D_z"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$3) w boxxy palette notitle 

set ylabel "y position (nm)"
set xlabel "x position (nm)"
set colorbox
set title "top layer, D_x"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($8/$4) w boxxy palette notitle  

unset colorbox 
unset ylabel 
unset xlabel 
set title "top layer, D_y"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$4) w boxxy palette notitle  

set title "top layer, D_z"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$4) w boxxy palette notitle  

unset multiplot 

file = "config_energy-1.1-3NN-DMI"
set output sprintf("%s.png", file)
set multiplot layout 2,3

set colorbox
set title "bottom layer, D_x"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

unset colorbox 
unset ylabel 
unset xlabel 

set title "bottom layer, D_y"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$3) w boxxy palette notitle 

set title "bottom layer, D_z"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$3) w boxxy palette notitle 

set ylabel "y position (nm)"
set xlabel "x position (nm)"
set colorbox
set title "top layer, D_x"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($8/$4) w boxxy palette notitle  

unset colorbox 
unset ylabel 
unset xlabel 
set title "top layer, D_y"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$4) w boxxy palette notitle  

set title "top layer, D_z"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$4) w boxxy palette notitle  

unset multiplot 

file = "config_energy-1.41-3NN-DMI"
set output sprintf("%s.png", file)
set multiplot layout 2,3

set colorbox
set title "bottom layer, D_x"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

unset colorbox 
unset ylabel 
unset xlabel 

set title "bottom layer, D_y"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$3) w boxxy palette notitle 

set title "bottom layer, D_z"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$3) w boxxy palette notitle 

set ylabel "y position (nm)"
set xlabel "x position (nm)"
set colorbox
set title "top layer, D_x"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($8/$4) w boxxy palette notitle  

unset colorbox 
unset ylabel 
unset xlabel 
set title "top layer, D_y"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$4) w boxxy palette notitle  

set title "top layer, D_z"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$4) w boxxy palette notitle  

unset multiplot 

file = "config_energy-1.41-2NN-DMI"
set output sprintf("%s.png", file)
set multiplot layout 2,3

set colorbox
set title "bottom layer, D_x"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

unset colorbox 
unset ylabel 
unset xlabel 

set title "bottom layer, D_y"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$3) w boxxy palette notitle 

set title "bottom layer, D_z"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$3) w boxxy palette notitle 

set ylabel "y position (nm)"
set xlabel "x position (nm)"
set colorbox
set title "top layer, D_x"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($8/$4) w boxxy palette notitle  

unset colorbox 
unset ylabel 
unset xlabel 
set title "top layer, D_y"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$4) w boxxy palette notitle  

set title "top layer, D_z"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$4) w boxxy palette notitle  

unset multiplot 

file = "config_energy-1.41-2NN-inter"
set output sprintf("%s.png", file)
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 

file = "config_energy-1.1-2NN-inter"
set output sprintf("%s.png", file)
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 


set xrange [2:198]
set yrange [2:198]
file = "config_energy-0.5-2NN-inter"
set output sprintf("%s.png", file)
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 

set cbrange [-1:1]
set xrange [2:198]
set yrange [2:198]
file = "config_energy-0.0-2NN-inter"
set output sprintf("%s.png", file)
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 

set cbrange [-0.1:0.1]
set xrange [2:198]
set yrange [2:198]
file = "config_energy-0.0-3NN-inter"
set output sprintf("%s.png", file)
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6) w boxxy palette notitle  

unset multiplot 

set cbrange [-0.1:0.1]
set xrange [2:198]
set yrange [2:198]
file = "config_energy-0.5-3NN-inter"
set output sprintf("%s.png", file)
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 

set xrange [2:98]
set yrange [2:98]
set cbrange [-0.1:0.1]
file = "config_energy-1.1-3NN-inter"
set output sprintf("%s.png", file)
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 

set cbrange [-0.1:0.1]
file = "config_energy-1.41-3NN-inter"
set output sprintf("%s.png", file)
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 


set cbrange [0.5:1.1]

set output "config_energy-1.41-1NN.png" 
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot "config_energy-1.41-1NN.txt" u 1:2:(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot "config_energy-1.41-1NN.txt" u 1:2:(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 

set cbrange [0.5:1.1]

set output "config_energy-1.1-1NN.png" 
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot "config_energy-1.1-1NN.txt" u 1:2:(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot "config_energy-1.1-1NN.txt" u 1:2:(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 

set output "config_energy-0.5-1NN.png" 
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot "config_energy-0.5-1NN.txt" u 1:2:(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot "config_energy-0.5-1NN.txt" u 1:2:(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 

set cbrange [0.54:1.0]
unset colorbox 
set output "config_energy-0.5-2NN.png" 
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot "config_energy-0.5-2NN.txt" u 1:2:(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot "config_energy-0.5-2NN.txt" u 1:2:(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 

set output "config_energy-1.1-2NN.png" 
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot "config_energy-1.1-2NN.txt" u 1:2:(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot "config_energy-1.1-2NN.txt" u 1:2:(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 

set cbrange [0.32:0.76]
set output "config_energy-1.41-3NN.png" 
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot "config_energy-1.41-3NN.txt" u 1:2:(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot "config_energy-1.41-3NN.txt" u 1:2:(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 

set output "config_energy-1.1-3NN.png" 
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot "config_energy-1.1-3NN.txt" u 1:2:(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot "config_energy-1.1-3NN.txt" u 1:2:(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 


set output "config_energy-0.5-3NN.png" 
set multiplot layout 1,2

set title "bottom layer"
set size 0.5,0.5
plot "config_energy-0.5-3NN.txt" u 1:2:(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "top layer"
set size 0.5,0.5
plot "config_energy-0.5-3NN.txt" u 1:2:(0.693):(0.6002):($6/$4) w boxxy palette notitle  

unset multiplot 


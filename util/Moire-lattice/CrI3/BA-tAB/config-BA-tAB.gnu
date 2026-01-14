

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


set term pngcairo size 800,800

set ytics 50 out 
set xtics 50 out 
set mytics 5 
set mxtics 5

set palette defined (  0 'blue', 1 'white',  2 'red')

set colorbox

set ytics
set xtics 
set format x ""
set format y ""

set xrange [2:198]
set yrange [2:198]


file = "config_energy-1.41-3NN"
set output sprintf("%s.png", file)
set multiplot layout 4,4

set  style fill solid noborder
set colorbox

set title "Layer 1, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($4/$3) w boxxy palette notitle 

set title "Layer 1, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "Layer 1, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$3) w boxxy palette notitle 

set title "Layer 1, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

set title "Layer 2, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$8) w boxxy palette notitle 

set title "Layer 2, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$8) w boxxy palette notitle 

set title "Layer 2, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$8) w boxxy palette notitle 

set title "Layer 2, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$8) w boxxy palette notitle 

set title "Layer 3, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($14/$13) w boxxy palette notitle 

set title "Layer 3, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($15/$13) w boxxy palette notitle 

set title "Layer 3, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($16/$13) w boxxy palette notitle 

set title "Layer 3, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($17/$13) w boxxy palette notitle 

set title "Layer 4, J"
plot file.".txt" u ($1*0.693):($2*0.600):(0.693):(0.600):($19/$18) w boxxy palette notitle 

set title "Layer 4, D_x"
plot file.".txt" u ($1*0.693):($2*0.600):(0.693):(0.6002):($20/$18) w boxxy palette notitle 

set title "Layer 4, D_y"
plot file.".txt" u ($1*0.693):($2*0.600):(0.693):(0.600):($21/$18) w boxxy palette notitle 

set title "Layer 4, D_z"
plot file.".txt" u ($1*0.693):($2*0.600):(0.693):(0.600):($22/$18) w boxxy palette notitle 

unset multiplot 

file = "config_energy-1.41-2NN"
set output sprintf("%s.png", file)
set multiplot layout 4,4

set colorbox

set title "Layer 1, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($4/$3) w boxxy palette notitle 

set title "Layer 1, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "Layer 1, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$3) w boxxy palette notitle 

set title "Layer 1, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

set title "Layer 2, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$8) w boxxy palette notitle 

set title "Layer 2, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$8) w boxxy palette notitle 

set title "Layer 2, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$8) w boxxy palette notitle 

set title "Layer 2, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$8) w boxxy palette notitle 

set title "Layer 3, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($14/$13) w boxxy palette notitle 

set title "Layer 3, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($15/$13) w boxxy palette notitle 

set title "Layer 3, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($16/$13) w boxxy palette notitle 

set title "Layer 3, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($17/$13) w boxxy palette notitle 

set title "Layer 4, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($19/$18) w boxxy palette notitle 

set title "Layer 4, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($20/$18) w boxxy palette notitle 

set title "Layer 4, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($21/$18) w boxxy palette notitle 

set title "Layer 4, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($22/$18) w boxxy palette notitle 

unset multiplot 


file = "config_energy-1.1-3NN"
set output sprintf("%s.png", file)
set multiplot layout 4,4

set colorbox

set title "Layer 1, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($4/$3) w boxxy palette notitle 

set title "Layer 1, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "Layer 1, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$3) w boxxy palette notitle 

set title "Layer 1, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

set title "Layer 2, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$8) w boxxy palette notitle 

set title "Layer 2, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$8) w boxxy palette notitle 

set title "Layer 2, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$8) w boxxy palette notitle 

set title "Layer 2, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$8) w boxxy palette notitle 

set title "Layer 3, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($14/$13) w boxxy palette notitle 

set title "Layer 3, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($15/$13) w boxxy palette notitle 

set title "Layer 3, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($16/$13) w boxxy palette notitle 

set title "Layer 3, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($17/$13) w boxxy palette notitle 

set title "Layer 4, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($19/$18) w boxxy palette notitle 

set title "Layer 4, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($20/$18) w boxxy palette notitle 

set title "Layer 4, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($21/$18) w boxxy palette notitle 

set title "Layer 4, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($22/$18) w boxxy palette notitle 

unset multiplot 

file = "config_energy-1.1-2NN"
set output sprintf("%s.png", file)
set multiplot layout 4,4

set colorbox

set title "Layer 1, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($4/$3) w boxxy palette notitle 

set title "Layer 1, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "Layer 1, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$3) w boxxy palette notitle 

set title "Layer 1, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

set title "Layer 2, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$8) w boxxy palette notitle 

set title "Layer 2, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$8) w boxxy palette notitle 

set title "Layer 2, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$8) w boxxy palette notitle 

set title "Layer 2, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$8) w boxxy palette notitle 

set title "Layer 3, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($14/$13) w boxxy palette notitle 

set title "Layer 3, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($15/$13) w boxxy palette notitle 

set title "Layer 3, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($16/$13) w boxxy palette notitle 

set title "Layer 3, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($17/$13) w boxxy palette notitle 

set title "Layer 4, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($19/$18) w boxxy palette notitle 

set title "Layer 4, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($20/$18) w boxxy palette notitle 

set title "Layer 4, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($21/$18) w boxxy palette notitle 

set title "Layer 4, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($22/$18) w boxxy palette notitle 

unset multiplot 

file = "config_energy-0.5-3NN"
set output sprintf("%s.png", file)
set multiplot layout 4,4

set colorbox

set title "Layer 1, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($4/$3) w boxxy palette notitle 

set title "Layer 1, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "Layer 1, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$3) w boxxy palette notitle 

set title "Layer 1, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

set title "Layer 2, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$8) w boxxy palette notitle 

set title "Layer 2, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$8) w boxxy palette notitle 

set title "Layer 2, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$8) w boxxy palette notitle 

set title "Layer 2, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$8) w boxxy palette notitle 

set title "Layer 3, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($14/$13) w boxxy palette notitle 

set title "Layer 3, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($15/$13) w boxxy palette notitle 

set title "Layer 3, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($16/$13) w boxxy palette notitle 

set title "Layer 3, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($17/$13) w boxxy palette notitle 

set title "Layer 4, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($19/$18) w boxxy palette notitle 

set title "Layer 4, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($20/$18) w boxxy palette notitle 

set title "Layer 4, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($21/$18) w boxxy palette notitle 

set title "Layer 4, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($22/$18) w boxxy palette notitle 

unset multiplot 

file = "config_energy-0.5-2NN"
set output sprintf("%s.png", file)
set multiplot layout 4,4

set colorbox

set title "Layer 1, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($4/$3) w boxxy palette notitle 

set title "Layer 1, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "Layer 1, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$3) w boxxy palette notitle 

set title "Layer 1, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

set title "Layer 2, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$8) w boxxy palette notitle 

set title "Layer 2, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$8) w boxxy palette notitle 

set title "Layer 2, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$8) w boxxy palette notitle 

set title "Layer 2, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$8) w boxxy palette notitle 

set title "Layer 3, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($14/$13) w boxxy palette notitle 

set title "Layer 3, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($15/$13) w boxxy palette notitle 

set title "Layer 3, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($16/$13) w boxxy palette notitle 

set title "Layer 3, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($17/$13) w boxxy palette notitle 

set title "Layer 4, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($19/$18) w boxxy palette notitle 

set title "Layer 4, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($20/$18) w boxxy palette notitle 

set title "Layer 4, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($21/$18) w boxxy palette notitle 

set title "Layer 4, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($22/$18) w boxxy palette notitle 

unset multiplot 


file = "config_energy-0.0-3NN"
set output sprintf("%s.png", file)
set multiplot layout 4,4

set colorbox

set title "Layer 1, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($4/$3) w boxxy palette notitle 

set title "Layer 1, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "Layer 1, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$3) w boxxy palette notitle 

set title "Layer 1, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

set title "Layer 2, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$8) w boxxy palette notitle 

set title "Layer 2, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$8) w boxxy palette notitle 

set title "Layer 2, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$8) w boxxy palette notitle 

set title "Layer 2, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$8) w boxxy palette notitle 

set title "Layer 3, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($14/$13) w boxxy palette notitle 

set title "Layer 3, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($15/$13) w boxxy palette notitle 

set title "Layer 3, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($16/$13) w boxxy palette notitle 

set title "Layer 3, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($17/$13) w boxxy palette notitle 

set title "Layer 4, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($19/$18) w boxxy palette notitle 

set title "Layer 4, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($20/$18) w boxxy palette notitle 

set title "Layer 4, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($21/$18) w boxxy palette notitle 

set title "Layer 4, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($22/$18) w boxxy palette notitle 

unset multiplot 

file = "config_energy-0.0-2NN"
set output sprintf("%s.png", file)
set multiplot layout 4,4

set colorbox

set title "Layer 1, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($4/$3) w boxxy palette notitle 

set title "Layer 1, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($5/$3) w boxxy palette notitle 

set title "Layer 1, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($6/$3) w boxxy palette notitle 

set title "Layer 1, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($7/$3) w boxxy palette notitle 

set title "Layer 2, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($9/$8) w boxxy palette notitle 

set title "Layer 2, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($10/$8) w boxxy palette notitle 

set title "Layer 2, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($11/$8) w boxxy palette notitle 

set title "Layer 2, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($12/$8) w boxxy palette notitle 

set title "Layer 3, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($14/$13) w boxxy palette notitle 

set title "Layer 3, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($15/$13) w boxxy palette notitle 

set title "Layer 3, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($16/$13) w boxxy palette notitle 

set title "Layer 3, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($17/$13) w boxxy palette notitle 

set title "Layer 4, J"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($19/$18) w boxxy palette notitle 

set title "Layer 4, D_x"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($20/$18) w boxxy palette notitle 

set title "Layer 4, D_y"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($21/$18) w boxxy palette notitle 

set title "Layer 4, D_z"
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693):(0.6002):($22/$18) w boxxy palette notitle 

unset multiplot 


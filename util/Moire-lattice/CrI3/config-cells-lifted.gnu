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


set ylabel "y position (nm)"
set xlabel "x position (nm)"

set ytics 50 in scale 2 mirror 
set xtics 50 in scale 2 mirror 
set mytics 5 
set mxtics 5 

chck(S,s, L, l, x) = (S == s && L == l) ? (x):(-100.0)
set colorbox
set cblabel "meV/mu_B^2"
set style fill solid noborder 
normalise = 2.98*2.98*1.602176634e-22

set palette defined (-0.15 '#81B1CB', 0.0 '#f5f0f7', 0.05 '#7A5FA9',  0.15'#614199')


do for [angle in "0.5"] {

file = sprintf("config-energy-cells-%s-AFMinter-DMIsub", angle) #, J_index, DMI_index, J_twist, J_prist, J_intra
print(file.".txt")
file = "config_energy_cells.txt"

set xrange [5:95]
set yrange [5:95]
set zrange [-11:11]
unset ztics
set auto cb 

set format x "" 
set format y ""

set ylabel ""
set xlabel ""

cell_x = 0.693*2
cell_y = 0.6002*2


set term pngcairo size 1200,1200
set output sprintf("%s.png", file)
set multiplot 
set pm3d at s interpolate 10,1
#set hidden3d 
set zeroaxis lw 0 

set view 60, 330
unset surface 
set contour base 
#unset cornerpoles 


set cbrange [0.06:-0.06]
set xyplane at -10
splot file u ($1*cell_x):($2*cell_y):(-10):(($27+$43)/($19+$35)/normalise) w pm3d lc palette notitle

set cbrange [0.12:0.07]
set xyplane at 0
splot file u ($1*cell_x):($2*cell_y):(0):(($26+$30+$42+$46+$25+$29+$41+$45+$24+$28+$40+$44)/($19+$35)/normalise/2.0) w pm3d lc palette notitle

#set cbrange [0.12:0.07]
#set xyplane at 10
#splot file u ($1*cell_x):($2*cell_y):(10):(($27+$43+($26+$30+$42+$46+$25+$29+$41+$45+$24+$28+$40+$44)/2.0)/($19+$35)/normalise) w pm3d lc palette notitle


}
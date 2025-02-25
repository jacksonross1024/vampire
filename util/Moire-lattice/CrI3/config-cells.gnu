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


set term pngcairo size 1600,2000


set ylabel "y position (nm)"
set xlabel "x position (nm)"

set ytics 50 out scale 2
set xtics 50 out scale 2
set mytics 5 
set mxtics 5 

set palette defined ( 0 'blue', 1 'white',  2 'red')

chck(S,s, L, l, x) = (S == s && L == l) ? (x):(-100.0)
set colorbox
set cblabel "meV/mu_B^2"
set style fill solid noborder 
normalise = 2.98*2.98

set xrange [5:195]
set yrange [5:195]

set palette defined (-0.15 '#81B1CB', 0.0 '#f5f0f7', 0.05 '#7A5FA9',  0.15'#614199')

do for [angle in "0.6"] {
do for [J_index in "0.5"] {
do for [DMI_index in "2"] {
do for [J_twist in "1.2"] {
do for [J_prist in "1.2"] {
do for [J_intra in "1.0"] {
file = sprintf("config_energy_cells") #, angle, J_index, DMI_index, J_twist, J_prist, J_intra
print(file.".txt")

cell_x = 0.693*2.0
cell_y = 0.6002*2.0

set view 75,45,1,1.
set pm3d at s interpolate 2,2
set output sprintf("%s-per-atom-12123434-lifted-Jinter.png", file)

set multiplot 
set xyplane 0
set zrange [0:14]
set format z ""

set xrange [5:195]
set yrange [5:195]

set format x "" 
set format y ""

set ylabel ""
set xlabel ""


unset colorbox 
set cbrange [0.4:-1]
set xyplane at 0
splot file.".txt" u ($1*cell_x):($2*cell_y):(0):(($13+$20)/($10+$17)) w pm3d at s notitle

set xyplane at 5
set cbrange [0.15:-0.15]
splot file.".txt" u ($1*cell_x):($2*cell_y):(5):(($14+$21)/($10+$17)) w pm3d at s notitle

set xyplane at 10
set cbrange [0.15:-0.15]
splot file.".txt"  u ($1*cell_x):($2*cell_y):(10):(($15+$22)/($10+$17)) w pm3d at s notitle


set xyplane at 14
splot file.".txt"  u ($1*cell_x):($2*cell_y):(14):(($16+$23)/($10+$17)) w pm3d at s notitle


unset multiplot


set output sprintf("%s-per-atom-12123434-lifted-Jintra-Jinter.png", file)

set multiplot 
set xyplane 0
set zrange [0:14]
set format z ""

set xrange [5:195]
set yrange [5:195]

set format x "" 
set format y ""

set ylabel ""
set xlabel ""


unset colorbox 
set cbrange [0.4:-1]
set xyplane at 0
splot file.".txt" u ($1*cell_x):($2*cell_y):(0):(($13+$20)/($10+$17)) w pm3d at s notitle

set xyplane at 5
set cbrange [18.3:17.4]
splot file.".txt" u ($1*cell_x):($2*cell_y):(5):(($12+$19)/($10+$17)) w pm3d at s notitle


unset multiplot

set output sprintf("%s-per-atom-2.png", file)
set multiplot layout 6,5


layer = 1234
set title sprintf("layer: %f; J_{intra}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):(($5+$12+$19+$26)/($3+$10+$17+$24)) w boxxy lc palette notitle

set title sprintf("layer: %f; J_{inter}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):(($6+$13+$20+$27)/($3+$10+$17+$24)) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{x}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):(($7+$14+$21+$28)/($3+$10+$17+$24)) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{y}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):(($8+$15+$22+$29)/($3+$10+$17+$24)) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{y}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):(($9+$16+$23+$30)/($3+$10+$17+$24)) w boxxy lc palette notitle


layer = 1
set title sprintf("layer: %f; J_{intra}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($5/$3) w boxxy lc palette notitle

set title sprintf("layer: %f; J_{inter}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($6/$3) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{x}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($7/$3) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{y}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($8/$3) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{y}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($9/$3) w boxxy lc palette notitle


layer = 2
set title sprintf("layer: %f; J_{intra}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($12/$10) w boxxy lc palette notitle

set title sprintf("layer: %f; J_{inter}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($13/$10) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{x}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($14/$10) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{y}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($15/$10) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{y}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($16/$10) w boxxy lc palette notitle


layer = 23
set title sprintf("layer: %f; J_{intra}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):(($12+$19)/($10+$17)) w boxxy lc palette notitle

set title sprintf("layer: %f; J_{inter}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):(($13+$20)/($10+$17)) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{x}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):(($14+$21)/($10+$17)) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{y}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):(($15+$22)/($10+$17)) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{y}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):(($16+$23)/($10+$17)) w boxxy lc palette notitle



layer = 3
set title sprintf("layer: %f; J_{intra}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($19/$17) w boxxy lc palette notitle

set title sprintf("layer: %f; J_{inter}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($20/$17) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{x}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($21/$17) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{y}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($22/$17) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{y}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($23/$17) w boxxy lc palette notitle


layer = 4
set title sprintf("layer: %f; J_{intra}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($26/$24) w boxxy lc palette notitle

set title sprintf("layer: %f; J_{inter}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($27/$24) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{x}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($28/$24) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{y}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($29/$24) w boxxy lc palette notitle

set title sprintf("layer: %f; D_{y}", layer)
plot file.".txt" u ($1*cell_x):($2*cell_y):(cell_x/2.0):(cell_y/2.0):($30/$24) w boxxy lc palette notitle

unset multiplot 


}
}
}
}
}
}

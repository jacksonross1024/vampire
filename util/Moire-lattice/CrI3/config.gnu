
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


set term pngcairo size 1600,1600


set ylabel "y position (nm)"
set xlabel "x position (nm)"

set ytics 100 out 
set xtics 100 out 
set mytics 50 
set mxtics 50

set palette defined ( 0 'blue', 1 'white',  2 'red')

chck(S,s, L, l, x) = (S == s && L == l) ? (x):(-100.0)
set colorbox
set cblabel "meV/mu_B^2"
set style fill solid noborder 
normalise = 2.98*2.98

set xrange [5:95]
set yrange [5:95]

set palette defined (-1.5 '#81B1CB', 0.0 '#f5f0f7', 0.5 '#7A5FA9',  1.5 '#614199')

do for [angle in "0.6"] {
do for [J_index in "0.5"] {
do for [DMI_index in "2"] {
do for [J_twist in "1.2"] {
do for [J_prist in "1.2"] {
do for [J_intra in "1.0"] {
file = sprintf("config_energy_atomic") #, angle, J_index, DMI_index, J_twist, J_prist, J_intra) 
print(file.".txt")



#set cbrange [0.6:-0.15]
set output sprintf("%s-J-net-twist-dmi-DMIsub.png", file)
set multiplot layout 4,4

layer = 0
species = 1
set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($8+$13+$18+$38+$43+$48)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($9+$14+$19+$39+$44+$49)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($10+$15+$20+$40+$45+$50)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($9+$14+$19+$39+$44+$49,$8+$13+$18+$38+$43+$48)) w boxxy palette notitle


species = 2
set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($8+$13+$18+$38+$43+$48)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($9+$14+$19+$39+$44+$49)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($10+$15+$20+$40+$45+$50)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($9+$14+$19+$39+$44+$49,$8+$13+$18+$38+$43+$48)) w boxxy palette notitle


#set cbrange [-0.15:0.6]
species = 3
set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($8+$13+$18+$38+$43+$48)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($9+$14+$19+$39+$44+$49)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($10+$15+$20+$40+$45+$50)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($9+$14+$19+$39+$44+$49,$8+$13+$18+$38+$43+$48)) w boxxy palette notitle


species = 4
set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($8+$13+$18+$38+$43+$48)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($9+$14+$19+$39+$44+$49)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($10+$15+$20+$40+$45+$50)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( ($8+$13+$18+$38+$43+$48) > 0.0 ? (0.3185*atan2( ($9+$14+$19+$39+$44+$49) ,($8+$13+$18+$38+$43+$48))) : (($9+$14+$19+$39+$44+$49) >= 0.0 ? (1+0.3185*atan2( ($9+$14+$19+$39+$44+$49) ,($8+$13+$18+$38+$43+$48))) : (-1+0.3185*atan2( ($9+$14+$19+$39+$44+$49) ,($8+$13+$18+$38+$43+$48)))  )) w boxxy palette notitle


unset multiplot 

#set cbrange [-0.08:-0.3]
set output sprintf("%s-J-intra-pristine-1nn-DMIsub.png", file)
set multiplot layout 4,5

layer = 0
species = 1
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($7/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($8/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($9/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($10/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($9,$8)) w boxxy palette notitle


layer = 0
species = 2
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($7/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($8/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($9/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($10/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($9,$8)) w boxxy palette notitle


layer = 0
species = 3
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($7/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($8/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($9/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($10/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($9,$8)) w boxxy palette notitle


layer = 0
species = 4
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($7/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($8/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($9/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($10/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($9,$8)) w boxxy palette notitle


unset multiplot 


set output sprintf("%s-J-intra-twist-1nn-DMIsub.png", file)
set multiplot layout 4,5

layer = 1
species = 1
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($7/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($8/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($9/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($10/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($9,$8)) w boxxy palette notitle


layer = 1
species = 2
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($7/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($8/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($9/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($10/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($9,$8)) w boxxy palette notitle



layer = 1
species = 3
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($7/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($8/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($9/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($10/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($9,$8)) w boxxy palette notitle


layer = 1
species = 4
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($7/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($8/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($9/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($10/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($9,$8)) w boxxy palette notitle


unset multiplot 


set output sprintf("%s-J-intra-pristine-2nn-DMIsub.png", file)
set multiplot layout 4,5

layer = 0
species = 1
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($12/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($13/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($14/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($15/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($14,$13)) w boxxy palette notitle


layer = 0
species = 2
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($12/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($13/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($14/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($15/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($14,$13)) w boxxy palette notitle


layer = 0
species = 3
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($12/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($13/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($14/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($15/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($14,$13)) w boxxy palette notitle


layer = 0
species = 4
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($12/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($13/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($14/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($15/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($14,$13)) w boxxy palette notitle


unset multiplot 


set output sprintf("%s-J-intra-twist-2nn-DMIsub.png", file)
set multiplot layout 4,5

layer = 1
species = 1
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($12/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($13/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($14/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($15/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($14,$13)) w boxxy palette notitle

layer = 1
species = 2
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($12/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($13/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($14/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($15/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($14,$13)) w boxxy palette notitle


layer = 1
species = 3
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($12/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($13/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($14/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($15/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($14,$13)) w boxxy palette notitle


layer = 1
species = 4
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($12/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($13/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($14/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($15/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($14,$13)) w boxxy palette notitle


unset multiplot 


set output sprintf("%s-J-intra-pristine-3nn-DMIsub.png", file)
set multiplot layout 4,5

layer = 0
species = 1
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($17/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($18/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($19/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($20/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($19,$18)) w boxxy palette notitle


layer = 0
species = 2
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($17/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($18/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($19/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($20/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($19,$18)) w boxxy palette notitle


layer = 0
species = 3
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($17/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($18/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($19/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($20/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($19,$18)) w boxxy palette notitle


layer = 0
species = 4
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($17/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($18/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($19/normalise) w boxxy palette notitle

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):($20/normalise) w boxxy palette notitle

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($19,$18)) w boxxy palette notitle


unset multiplot 


set output sprintf("%s-J-intra-twist-3NN-DMIsub.png", file)
set multiplot layout 4,5

layer = 1
species = 1
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($17)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($18)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($19)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($20)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($19,$18)) w boxxy palette notitle


layer = 1
species = 2
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($17)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($18)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($19)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($20)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($19,$18)) w boxxy palette notitle


layer = 1
species = 3
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($17)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($18)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($19)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($20)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($19,$18)) w boxxy palette notitle


layer = 1
species = 4
set title sprintf("species %.f, layer %.f, J_{intra}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($17)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($18)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($19)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($20)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($19,$18)) w boxxy palette notitle

unset multiplot 



set output sprintf("%s-J-inter-pristine-DMIsub.png", file)
set multiplot layout 4,1

layer = 0
species = 1
set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($22+$27+$32)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_x", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($23+$28+$33)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_y", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($24+$29+$34)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_z", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($25+$30+$35)/normalise) w boxxy palette notitle 

layer = 0
species = 2
set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($22+$27+$32)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($23+$28+$33)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_y", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($24+$29+$34)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_z", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($25+$30+$35)/normalise) w boxxy palette notitle 


layer = 0
species = 3
set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($22+$27+$32)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_x", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($23+$28+$33)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_y", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($24+$29+$34)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_z", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($25+$30+$35)/normalise) w boxxy palette notitle 


layer = 0
species = 4
set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($22+$27+$32)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($23+$28+$33)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_y", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($24+$29+$34)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_z", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($25+$30+$35)/normalise) w boxxy palette notitle 


unset multiplot 


# set output sprintf("%s-J-inter-twist-DMIsub.png", file)
# set multiplot layout 4,4

# layer = 1
# species = 1
# set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($22+$27+$32)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_x", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($23+$28+$33)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_y", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($24+$29+$34)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_z", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($25+$30+$35)/normalise) w boxxy palette notitle 


# species = 2
# set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($22+$27+$32)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_x", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($23+$28+$33)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_y", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($24+$29+$34)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_z", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($25+$30+$35)/normalise) w boxxy palette notitle 


# species = 3
# set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($22+$27+$32)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_x", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($23+$28+$33)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_y", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($24+$29+$34)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_z", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($25+$30+$35)/normalise) w boxxy palette notitle 


# species = 4
# set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($22+$27+$32)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_x", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($23+$28+$33)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_y", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($24+$29+$34)/normalise) w boxxy palette notitle 

# set title sprintf("layer %.f, D_z", layer)
# plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($25+$30+$35)/normalise) w boxxy palette notitle 


# unset multiplot 

#set cbrange [0.6:-0.15]
set output sprintf("%s-J-inter-twist-2-DMIsub.png", file)
set multiplot layout 4,5

layer = 1
species = 1
set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($37+$42+$47)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($38+$43+$48)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($39+$44+$49)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($40+$45+$50)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($39+$44+$49,$38+$43+$48)) w boxxy palette notitle


species = 2
set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($37+$42+$47)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($38+$43+$48)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($39+$44+$49)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($40+$45+$50)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($39+$44+$49,$38+$43+$48)) w boxxy palette notitle


#set cbrange [-0.15:0.6]
species = 3
set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($37+$42+$47)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($38+$43+$48)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($39+$44+$49)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($40+$45+$50)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($39+$44+$49,$38+$43+$48)) w boxxy palette notitle


species = 4
set title sprintf("species %.f, layer %.f, J_{inter}", species, layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($37+$42+$47)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($38+$43+$48)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($39+$44+$49)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):(($40+$45+$50)/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, theta_{Dxy}", layer)
plot file.".txt" u (chck($2,species, $3, layer, $4*0.1)):($5*0.1):(0.693/2.0):(0.6002/2.0):( 0.3185*atan2($39+$44+$49,$38+$43+$48)) w boxxy palette notitle


unset multiplot 



}
}
}
}
}
}

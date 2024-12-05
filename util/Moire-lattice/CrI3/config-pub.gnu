
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

set ytics 50 out 
set xtics 50 out 
set mytics 5 
set mxtics 5

set palette defined ( 0 'blue', 1 'white',  2 'red')

chck(S,s,E,n) = (S == s) ? (E/n):(1/0)
set colorbox
set cblabel "meV/mu_B^2"
set style fill solid noborder 

set xrange [2:198]
set yrange [2:198]

do for [angle in "0.5 1.1 2.0"] {
do for [J_index in "0.2 0.25 0.3"] {
do for [DMI_index in "7"] {
do for [J_twist in "0.98 0.96 0.94 0.92"] {
#do for [J_twist in "0.95"] {

file = sprintf("config-energy-cells-%s-%s-%s-%s", angle, J_index, DMI_index, J_twist) 
print(file.".txt")
set output sprintf("%s.png", file)


set multiplot layout 4,4

normalise = 2.98*2.98

layer = 1
#set title sprintf("layer %.f, J_{inter}", layer)
#plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($5/normalise) w boxxy palette notitle 

#set title sprintf("layer %.f, D_x", layer)
#plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($6/normalise) w boxxy palette notitle 

#set title sprintf("layer %.f, D_y", layer)
#plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($7/normalise) w boxxy palette notitle 

#set title sprintf("layer %.f, D_z", layer)
#plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($8/normalise) w boxxy palette notitle 

layer = 2
set title sprintf("layer %.f, J", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($10/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($11/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($12/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($13/normalise) w boxxy palette notitle 

layer = 3
set title sprintf("layer %.f, J", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($15/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($16/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($17/normalise) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($18/normalise) w boxxy palette notitle 


set title "layer 2 + 3"
#set title sprintf("layer %s, J", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):(($10+$15)/normalise) w boxxy palette notitle 

#set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):(($11+$16)/normalise) w boxxy palette notitle 

#set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):(($12+$17)/normalise) w boxxy palette notitle 

#set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):(($13+$18)/normalise) w boxxy palette notitle 

set title "layer 1 + 2 + 3 + 4"
#set title sprintf("layer %s, J", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):(($5+$10+$15+$20)/normalise) w boxxy palette notitle 

#set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):(($6+$11+$16+$21)/normalise) w boxxy palette notitle 

#set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):(($7+$12+$17+$22)/normalise) w boxxy palette notitle 

#set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):(($8+$13+$18+$23)/normalise) w boxxy palette notitle 


#set title sprintf("layer %.f, J", layer)
#plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($20/normalise) w boxxy palette notitle 

#set title sprintf("layer %.f, D_x", layer)
#plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($21/normalise) w boxxy palette notitle 

#set title sprintf("layer %.f, D_y", layer)
#plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($22/normalise) w boxxy palette notitle 

#set title sprintf("layer %.f, D_z", layer)
#plot file.".txt" u ($1*0.693):($2*0.6002):(0.693/2.0):(0.6002/2.0):($23/normalise) w boxxy palette notitle 


unset multiplot 
}
}
}
}

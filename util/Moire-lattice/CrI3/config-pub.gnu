
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


set term pngcairo size 1600,1600 font "helvetica, 14"


unset ylabel 
unset xlabel 

unset ytics
unset xtics 

set cbrange [-0.15:0.6]


chck(S,s,E,n) = (S == s) ? (E/n):(1/0)
set colorbox user origin 0.8,0.25 size 0.025,0.3
set size 0.8,0.8
set style fill solid 0.95 noborder 

set xrange [2:198]
set yrange [2:198]

angle = '1.1'
J_index = '0.2'
DMI_index = '1.0'
J_twist = '1.0'
J_intra = '1.0'

set ticslevel 0  
set dgrid3d 30,30  
#set palette defined (0 "blue", 0.75 "white", 1.4 "red")  
set style lines 100 lt 5 lw 0.5  
#set pm3d hidden3d 100  
set grid  
set view 74,216  
unset key  

do for [angle in "0.5 1.1 2.0"] {
do for [J_index in "0.2 0.25 0.3 0.35"] {
do for [J_twist in "1.0 0.9 0.8"] {
do for [J_intra in "1.0"] {
do for [DMI_index in "2 4 8"] {

set palette defined (-0.15 '#81B1CB', 0.0 '#f5f0f7', 0.1 '#7A5FA9',  0.6'#614199')

file = sprintf("config-energy-cells-%s-%s-%s-%s-%s-Jinter", angle, J_index, DMI_index, J_twist, J_intra) 
print(file.".txt")
normalise = 2.98*2.98*2

set output sprintf("%s.png", file)
set cbrange [-0.15:0.6]
#set auto cb
set cbtics -0.15, 0.35,0.6 nomirror scale 1.0
set cbtics add  ("0" 0.0)
set cbtics add  ("0.6" 0.6)
set zrange [-0.15:0.6]

clr(c,e) = (c > 24 || c < 24) ? ( 1/0) : ( e)


cell = 1.0
len(x,y,z) = sqrt(x*x + y*y + z*z)

splot file.".txt" u (clr($9,$2*0.6)):($1*0.693):(($10+$15)/normalise) w pm3d notitle 
#show cbrange
}
}
}
}
}
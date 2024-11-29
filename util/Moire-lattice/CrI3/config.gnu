
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

set palette defined (  0 'blue', 1 'white',  2 'red')

chck(S,s,E,n) = (S == s) ? (E/n):(1/0)
set colorbox
set style fill solid noborder 

set xrange [2:198]
set yrange [2:198]

file = "config_energy"
set output sprintf("%s.png", file)
set multiplot layout 5,4



layer = 2
set title sprintf("layer %.f, J", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $4,$8)) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $5,$8))  w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $6,$8)) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $7,$8)) w boxxy palette notitle  

layer = 1
set title sprintf("layer %.f, J", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $4,$8)) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $5,$8))  w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $6,$8)) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $7,$8)) w boxxy palette notitle  

layer = 4
set title sprintf("layer %.f, J", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $4,$8)) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $5,$8))  w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $6,$8)) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $7,$8)) w boxxy palette notitle  

layer = 3
set title sprintf("layer %.f, J", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $4,$8)) w boxxy palette notitle 

set title sprintf("layer %.f, D_x", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $5,$8))  w boxxy palette notitle 

set title sprintf("layer %.f, D_y", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $6,$8)) w boxxy palette notitle 

set title sprintf("layer %.f, D_z", layer)
plot file.".txt" u ($2*0.1):($3*0.1):(0.693):(0.6002):(chck($1,layer, $7,$8)) w boxxy palette notitle  


unset multiplot 

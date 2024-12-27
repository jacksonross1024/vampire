

set palette defined (0 "blue", 1 "white", 2 "red")

set term pngcairo size 400, 300

set style fill solid noborder

set ylabel "Position (nm)"
set xlabel "Position (nm)"

set colorbox 


cell_x = 5
cell_y = 5
cx = 0.1
cy = 0.1

#set cbrange [-1:1]
set yrange [0:100]
set xrange [0:100]
file_index = 1

set output "corr-data.png"
plot "corr-dipole.txt" u ($1*cell_x):($2*cell_y):(cell_x/2):(cell_y/2):($3/1.6e19 - 1) w boxxy palette notitle 



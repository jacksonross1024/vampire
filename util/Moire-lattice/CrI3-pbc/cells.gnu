

set palette defined (0 "blue", 1 "white", 2 "red")

set term pngcairo 

set style fill solid noborder

set ylabel "Position (nm)"
set xlabel "Position (nm)"

set colorbox 
set cblabel "m_z (u_B/nm^2)"

cell_x = 2
cell_y = 2
cx = 0.2
cy = 0.2

set yrange [0:]
set xrange [0:]
file_index = 9
file = sprintf("cells-%08.f.txt", file_index)
print file 
set output sprintf("cells-%08.f-whole-app.png", file_index)
plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($26*$27/(cell_x*cell_y)) w boxxy palette notitle 

set output sprintf("cells-%08.f-1-app.png", file_index)
plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($6*$7*($8/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 

set output sprintf("cells-%08.f-2-app.png", file_index)
plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($11*$12*($13/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 

set output sprintf("cells-%08.f-3-app.png", file_index)
plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($16*$17*($18/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 

set output sprintf("cells-%08.f-4-app.png", file_index)
plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($21*$22*($23/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 



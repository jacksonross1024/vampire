



set palette defined (0 "blue", 1 "white", 2 "red")

set term pngcairo size 1200,2000

set style fill solid noborder

set ylabel "Position (nm)"
set xlabel "Position (nm)"

set colorbox 
set cblabel "m_z (u_B/nm^2)"
set cbrange [-1:1]
cell_x = 2
cell_y = 2
cx = 0.1
cy = 0.1

set yrange [0:60]
set xrange [0:50]
file_index = 20
folder = "xcross-2-AFM"
file = sprintf("cells-%08.f.txt", file_index)

set output sprintf("cells-%08.f-%s.png", file_index, folder)
set multiplot layout 5,3

plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($4*$7) w boxxy palette notitle 
plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($5*$7) w boxxy palette notitle 
plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($6*$7) w boxxy palette notitle 

plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($12*$9) w boxxy palette notitle 
plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($12*$10) w boxxy palette notitle 
plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($12*$11) w boxxy palette notitle 

plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($27*$24*$28*3.1/(cell_x*cell_y)) w boxxy palette notitle 
plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($27*$25*$28*3.1/(cell_x*cell_y)) w boxxy palette notitle 
plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($27*$26*$28*3.1/(cell_x*cell_y)) w boxxy palette notitle 

plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($14*$17) w boxxy palette notitle 
plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($15*$17) w boxxy palette notitle 
plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($16*$17) w boxxy palette notitle 

plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($19*$22) w boxxy palette notitle 
plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($20*$22) w boxxy palette notitle 
plot file u ($1*cx):(cy*$2):(cell_x/2.0):(cell_y/2.0):($21*$22) w boxxy palette notitle 


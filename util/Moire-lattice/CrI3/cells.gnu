


set palette defined (0 "blue", 1 "white", 2 "red")

set term pngcairo size 1600, 2500

set style fill solid noborder

set ylabel "Position (nm)"
set xlabel "Position (nm)"

set colorbox 


cell_x = .99
cell_y = .99
cx = 0.1
cy = 0.1

set yrange [0:500]
set xrange [0:500]
file_index = 1


#do for [temp_index in "5 15 25"]  {
#    do for [j_index in "0.2 0.225 0.25 0.275"]  {
j_index = 0.2
d_index = 10.0
do for [index=0:9]{
    file = sprintf("cells-%08.f.txt", index)
    set output sprintf("cells-5K-%.f-%.fJ-%.fDMI.png", index, j_index, d_index)
    print file 

    set multiplot layout 5,3

    set cblabel "m_x (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($19*$22*($23/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_y (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($20*$22*($23/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($21*$22*($23/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 

    #top
    set cblabel "m_x (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($14*$17*($18/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_y (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($15*$17*($18/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($16*$17*($18/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 
    
    #combined
    set cblabel "m_x (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($24*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_y (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($25*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($26*$27/(cell_x*cell_y)) w boxxy palette notitle 

    #bottom
    set cblabel "m_x (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($9*$12*($13/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_y (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($10*$12*($13/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($11*$12*($13/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle

    set cblabel "m_x (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($4*$7*($8/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_y (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($5*$7*($8/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($6*$7*($8/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle
    
    unset multiplot
}
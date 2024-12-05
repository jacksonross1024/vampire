


set palette defined (0 "blue", 1 "white", 2 "red")

set term pngcairo size 1600, 2000

set style fill solid noborder

set ylabel "Position (nm)"
set xlabel "Position (nm)"

set colorbox 



cx = 0.1
cy = 0.1

set yrange [0:500]
set xrange [0:500]
index = 1
d_index = 10.0
J_index = 0.2
angle = "0.2"
DMI_index = 4.0

do for [angle in "0.5 1.1 2.0"] {
do for [J_index in "0.2 0.25 0.3"] {
do for [DMI_index in "4 7 10"] {
do for [J_twist in "0.98 0.96 0.94 0.92"] {
    cell_x = 4.99
    cell_y = 4.99
    file = sprintf("dipole-cells-%s-5K-%s-%s-%s-AF.data", angle, J_index, DMI_index, J_twist)
    print file 
    set output sprintf("cells-%s-5K-dipole-%sJ-%sDMI-%s-AF.png", angle, J_index, DMI_index, J_twist)

    normalise = 1e-18
    set multiplot layout 4,3

    set cblabel "m_z (u_B/nm^2)"
    plot file u ($4*cx):(cy*$5):(cell_x):(cell_y):($6 == 19.62 ? $18*normalise : NaN) w boxxy palette notitle 
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($4*cx):(cy*$5):(cell_x):(cell_y):($6 == 19.62 ? $19*normalise : NaN) w boxxy palette notitle 
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($4*cx):(cy*$5):(cell_x):(cell_y):($6 == 19.62 ? $20*normalise : 1/0) w boxxy palette notitle 

    cell_x = 5
    cell_y = 5

    file = sprintf("cells-%s-5K-%s-%s-%s-AF.txt", angle, J_index, DMI_index, J_twist)
   # file = "cells-00000001.txt"
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($19*$22*($23/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_y (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($20*$22*($23/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($21*$22*($23/$28)*$27/(cell_x*cell_y)) w boxxy palette notitle 

    #combined
    set cblabel "m_x (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($24*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_y (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($25*$27/(cell_x*cell_y)) w boxxy palette notitle 
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):($26*$27/(cell_x*cell_y)) w boxxy palette notitle 

    w_1 = 0.41
    w_2 = 0.269
    w_3 = 0.186
    w_4 = 0.183
    #weighted
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):((w_1*$19*$22*$23+w_2*$14*$17*$18+w_3*$9*$12*$13+w_4*$4*$7*$8)*$27/(cell_x*cell_y*$28)) w boxxy palette notitle 
    set cblabel "m_y (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):((w_1*$20*$22*$23+w_2*$15*$17*$18+w_3*$10*$12*$13+w_4*$5*$7*$8)*$27/(cell_x*cell_y*$28)) w boxxy palette notitle 
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($1*cx):(cy*$2):(cell_x):(cell_y):((w_1*$21*$22*$23+w_2*$15*$17*$18+w_3*$11*$12*$13+w_4*$6*$7*$8)*$27/(cell_x*cell_y*$28)) w boxxy palette notitle 

unset multiplot 
}
}
}
}
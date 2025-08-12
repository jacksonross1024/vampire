set palette defined (0 "blue", 1 "white", 2 "red")

set term pngcairo size 1600, 2000

set style fill solid noborder

set ylabel "Position (nm)"
set xlabel "Position (nm)"

set colorbox 


cx = 0.1
cy = 0.1

set yrange [0:100]
set xrange [0:100]
index = 1
d_index = 10.0
J_index = 0.2
angle = "0.2"
DMI_index = 4.0

do for [angle in "1.1"] {
do for [J_index in "0.25"] {
do for [DMI_index in "7"] {
do for [J_twist in "1.0"] {
do for [J_intra in "1.0"] {
    cell_x = 2.5
    cell_y = 2.5
    file = sprintf("dipole-cells-test.data")
    print file 
    set output sprintf("dipole-cells-test.png")

    normalise = 1e-18
    set multiplot layout 4,3

    set cblabel "m_z (u_B/nm^2)"
    plot file u ($6 == 13.08 ? $4*cx : -200):(cy*$5):(cell_x):(cell_y):($20*normalise) w boxxy palette notitle 
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($6 == 19.62 ? $4*cx : -200):(cy*$5):(cell_x):(cell_y):($6 == 19.62 ? $20*normalise : NaN) w boxxy palette notitle 
    set cblabel "m_z (u_B/nm^2)"
    plot file u ($6 == 60 ? $4*cx : -200):(cy*$5):(cell_x):(cell_y):($6 == 60 ? $20*normalise : 1/0) w boxxy palette notitle 

    cell_x = 0.5
    cell_y = 0.5

    file = sprintf("cells-00000001.txt")
   # file = "cells-00000000.txt"
    #file = "cells-00000001.txt"
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
    plot file u ($3 == 0 ? ($1*cx) : -50):(cy*$2):(cell_x):(cell_y):($31*$32/(cell_x*cell_y)) w boxxy palette notitle 

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
}
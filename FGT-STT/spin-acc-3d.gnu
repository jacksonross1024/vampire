

set term wxt
set output "js_bottom.png" 

set xlabel "x position (nm)"
set ylabel "y position (nm)"

set mouse

unset ztics
set colorbox 
set auto cb 

set pm3d interpolate 2,2
set pm3d lighting primary 0.50 specular 0.2
set view 20,50

len(x,y,z) = sqrt(x*x + y*y + z*z)
splot "spin-acc/1-b" u ($3*0.1*3.99):($2*1):(len($13, $14, $15)) w pm3d notitle

set cbrange [0.997:1.00]
m_inf = 0.7293e7
splot "spin-acc/1-b" u ($3*0.1*3.99):($2*1):(len($7, $8, $9)/m_inf) w pm3d notitle


set pm3d spotlight rgb "blue" rot_x 90 Phong 0.1
set pm3d lighting primary 0.50 spec 0.2 spec2 0.2

set auto cb
splot "spin-acc/1-b" u ($3*0.1*3.99):($2*1):($26) w pm3d notitle
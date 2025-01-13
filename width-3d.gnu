#set global styles
# General Palette
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

# Reset for 3D/map plots
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


chck(sl_y, sl_z, mag, y, z) = (sl_y == y && sl_z == z) ? (mag) : (1/0)

set terminal pngcairo font "Helvetica, 14"

set output "torque.png"
set auto y
set auto x
s = 1e-4 
set ylabel "Torque (J)"
set xlabel "position (nm)"
plot "output" u ($1*s):18 w l title "T_x",\
"" u ($1*s):19 w l title "T_y",\
"" u ($1*s):20 w l title "T_z" 


p_2 = 250
w_2 = 20
m_y_2 = 1.0/sqrt(2.0)


#K = ARG1 + 0
K = 0.0
m_y_2 = m_y_2*(1.0-(K/1225.0)**2.0)**0.332


wdth_2(x, w, p,m ) = m*(sin(atan(exp(1.0*(x-p)/(w)))) - cos(atan(exp(1.0*(x-p)/(w))))) 

file = 1
file_d = sprintf("dw/dw-%.f.txt", file)

set fit quiet 
set fit errorvariables 

set auto x 

set print "dw-mc-dynamics.txt"
files = system("ls dw/dw-*.txt | sort --version-sort")
do for [file in files] {

do for [y=0:0] {
    do for [z=0:0] {
        
        fit wdth_2(x,w_2,p_2,m_y_2) file  u ($1*0.1):(chck($2,$3,$5,y,z)) via w_2, p_2, m_y_2
        
       print sprintf("%s %f %f %f %f %f %f %f %f", file, y, z,p_2,w_2,m_y_2,p_2_err,w_2_err,m_y_2_err)

    }
}
}


set output "position-dynamic.png"
set auto x 
set auto y 
set auto y2 
set y2tics out nomirror 
set ytics out nomirror 
set xtics out nomirror 
set ylabel "Position (nm)" tc "#20547C"
set xlabel "Time (ps)"
set y2label "Velocity (km/s)" tc "#C07F28"

x1 = 249.916804
x2 = 249.916803
dx = 0.0
d_x(x) = (x2 = x, dx = (x2-x1)/0.1, x1 = x2, dx) #nm/ps
set y2range [0:5]

plot "dw-mc-dynamics.txt" u ($0*0.1+2):4 w l ls 3 title "position",\
"" u ($0*0.1+2):(d_x($4)) axis x1y2 w l ls 6 title "Velocity",\


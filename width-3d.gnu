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

w_2 = 20.0
w_1 = 20.0
p_2 = 300
p_1 = 300

mag1 = sqrt(2.0)*0.5
chck(sl_y, sl_z, mag, y, z) = (sl_y != y && sl_z != z) ? (1/0) : (mag)

wdth_180_exp(x, w, p ) = cos(2.0*atan(exp((x-p)/w))) #k2|| <or> k4||
wdth_180_sinh(x, w, p ) = sin(-atan(sinh((x-p)/w))) #k2|| <and> k4||

wdth_90_s(x, w, p ) = sin(atan(exp(-(x-p)/w))-0.25*3.14) #k4||
wdth_90_c(x, w, p ) = cos(atan(exp(-3.14*0.5*(x-p)/w))-0.25*3.14) #k4||

l = 0.3328
p_1 = l*1791.0
p_2 = l*1791.0
p_1 = 300.0
p_2 = 300.0
set terminal pngcairo 

set output "torque.png"
set auto y
set auto x
s = 1e-4 
set ylabel "Torque (J)"
set xlabel "position (nm)"
plot "output" u ($1*s):40 w l title "T_x",\
"" u ($1*s):40 w l title "T_y",\
"" u ($1*s):40 w l title "T_z" 

p_1 = 501
w_1 = 86
p_2 = p_1
w_2 = w_1
m_y_2 = 1.0/sqrt(2.0)
p_1_t = 501
w_1_t = 86
p_2_t = p_1
w_2_t = w_1
m_y_t_2 = 1.0/sqrt(2.0)

K = ARG1 + 0

m_y_1 = m_y_2*(1.0-(K/1225.0)**2.0)**0.332
m_y_t_1 = m_y_1

wdth_2(x, w, p,m ) = m*(sin(atan(exp(1.0*(x-p)/(w)))) - cos(atan(exp(1.0*(x-p)/(w))))) 
wdth_tanh(x, w, p,m ) = m*tanh((x-p)/w)

file = 1
file_d = sprintf("dw/dw-%.f.txt", file)

set fit quiet 
#set print sprintf("%.0f-3d-data.txt", file)
set fit errorvariables 

set auto x 

set print "dw-mc-dynamics.txt"
files = system("ls dw/dw-*.txt | sort --version-sort")
do for [file in files] {

set print "dw-temp-data.txt"
do for [y=0:0] {
    do for [z=0:0] {
        fit wdth_2(x,w_1,p_1,m_y_1) file  u ($1*0.1):(chck($2,$3,-$5,y,z)) via w_1, p_1
        fit wdth_2(x,w_2,p_2,m_y_2) file  u ($1*0.1):(chck($2,$3,-$5,y,z)) via w_2, p_2, m_y_2
        fit wdth_tanh(x,w_1_t,p_1_t,m_y_t_1) file  u ($1*0.1):(chck($2,$3,-$5,y,z)) via w_1_t, p_1_t
        fit wdth_tanh(x,w_2_t,p_2_t,m_y_t_2) file  u ($1*0.1):(chck($2,$3,-$5,y,z)) via w_2_t, p_2_t, m_y_t_2
       print sprintf("%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", file, y, z,p_1,w_1,p_1_err,w_1_err, p_2, w_2, m_y_2, p_2_err, w_2_err, m_y_2_err,\
                                                                          p_1_t,2*w_1_t/3.14,p_1_t_err,w_1_t_err,p_2_t,w_2_t*2/3.14,m_y_t_2,p_2_t_err,w_2_t_err,m_y_t_2_err)
    }
}
stats "dw-temp-data.txt" u 5 name 'SG_1_'
stats "dw-temp-data.txt" u 9 name 'SG_2_'
stats "dw-temp-data.txt" u 15 name 'tan_1_'
stats "dw-temp-data.txt" u 19 name 'tan_2_'
set print "dw-mc-dynamics.txt" append
print SG_1_mean, SG_1_stddev, SG_2_mean, SG_2_stddev, tan_1_mean, tan_1_stddev, tan_2_mean, tan_2_stddev, p_2
}


set output "width-dynamic.png"
set auto x 
set ylabel "width (nm)" 
set xlabel "Time Step"


plot "dw-mc-dynamics.txt" u ($0*100.0):1 w l ls 1 title "m_e(0) SG",\
"" u ($0*100.0):3 w l ls 4 title "m_e(T) SG",\
"" u ($0*100.0):5 w l ls 5 title "m_e(0) tanh",\
"" u ($0*100.0):7 w l ls 7 title "m_e(T) tanh",\

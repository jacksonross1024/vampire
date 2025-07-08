
set term pngcairo size 800,600

set xlabel "Temperature (K)"
set ylabel "magnetisation (m_B)"

set ytics out  nomirror
set xtics out nomirror

set yrange [-0.1:2.75]
set ytics 0,0.5
set mytics 5

T_c_1 = 500
T_c_2 = 500
T_c_3 = 430
T_c_4 = 500
T_c_5 = 367
alpha = 1.0
beta = 0.323

alpha_1 = 1.45
beta_1 = 0.323
beta_2 = 0.323
m(T, a, b, T_c) = (T <= T_c) ? (1-(T/T_c)**a)**b : 0.0

#fit [0:T_c_1]  m(x,alpha,beta, T_c_1) "output-Li" u 1:2 via T_c_1


set output "Fe5GeTe2-mvsT.png"
set multiplot 

set xrange [0:600]
set key top right title "Ershadrad et al."
plot "Fe5GeTe2-Fe-mvsT.csv" u 1:10 w lp ls 1 title "Fe_5",\
"Fe5GeTe2-Fe-mvsT.csv" u 1:8 w lp ls 2 title "Fe_4",\
"Fe5GeTe2-Fe-mvsT.csv" u 1:6 w lp ls 3 title "Fe_3",\
"Fe5GeTe2-Fe-mvsT.csv" u 1:4 w lp ls 4 title "Fe_2",\
"Fe5GeTe2-Fe-mvsT.csv" u 1:2 w lp ls 5 title "Fe_1"

set key center right title "MC" 
plot "output" u 1:($5*0.109) w lp ls 1 dt " - " title "Fe_1",\
"" u 1:($9*2.288) w lp ls 2 dt " - " title "Fe_2",\
"" u 1:($13*2.055) w lp ls 3 dt " - " title "Fe_3",\
"" u 1:($17*1.458) w lp ls 4 dt " - " title "Fe_4",\
"" u 1:($21*2.5634) w lp ls 5 dt " - "  title "Fe_5"

unset  multiplot
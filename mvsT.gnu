
set term pngcairo size 800,600

set xlabel "Temperature (K)"
set ylabel "magnetisation (1/m_s)"

set ytics out  nomirror
set xtics out nomirror

set yrange [-0.1:1.1]
set ytics 0,0.25
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

fit [0:T_c_1]  m(x,alpha,beta, T_c_1) "output-Li" u 1:2 via T_c_1
fit [0:T_c_2] m(x,alpha,beta_1, T_c_2) "output-Li-rescaled" u 1:2 via T_c_2, alpha
fit  [0:T_c_5] m(x,alpha_1,beta_2, T_c_5) "Fe3GaTe2-exp-mvsT.csv" u 1:2 via T_c_5


set output "Fe3GaTe2-mvsT.png"
plot "Fe3GaTe2-exp-mvsT.csv" u 1:2 w p ls 1 title "Exp",\
m(x,alpha_1, beta_2, T_c_5) w l ls 1 title sprintf("T_c = %.0f, alpha = %.2f, beta = %.3f", T_c_5, alpha_1, beta_2),\
(0) w l dt " - " lc "grey" notitle,\
"output-lda-lr" u 1:2 w p ls 3 title "LDA, LR",\
"output-Li" u 1:2 w p ls 4 title "Li",\
"output-Li-rescaled" u 1:2 w p ls 5 title "Li rescaled",\
m(x,alpha, beta_1, T_c_2) w l ls 5 title sprintf("T_c = %.0f, alpha = %0.2f, beta = %.3f", T_c_2, alpha, beta_1)

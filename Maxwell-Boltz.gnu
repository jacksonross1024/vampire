#!/bin/bash

Ep = 300.0*652.0/6020.0
Ee = 300.0*252.0/6020.0
e  = 105.59
n  = 8788
w  = 2.0/8788.0

mb(x,t) = 0.15*sqrt(2.0/3.141592)*x*x*exp(-0.5*x*x/(t*t))/(t*t*t)
mb_e(x,t) = sqrt(2.0/3.141592)*x*x*exp(-1.0*x/t)/(t**1.5)
mb_TH(x,t) = (1.0/4.0)*sqrt(1*x/3.141592)*((1/t)**1.5)*exp(-x/t)
mb_S(x,t) = (1.0/120.0)*sqrt(1/(2.0*3.141592))*((1.0/t)**6.0)*(x**5.0)*exp(-1.0*x/t)

bin(x,w) = w*floor(x/w)


set term wxt 1
set title "Lattice parameter distribution comparison"
p "P_distrib" u ($1*2.0/n):(mb($1*w,Ep/e)) w l lw 2 title "MB",\
"" u ($1*2.0/n):(mb_e($1*w,Ep/e)) w l lw 2 title "MB_e",\
"" u ($1*2.0/n):(mb_TH($1*w,Ep/e)) w l lw 2 title "MB-1.5",\
"" u ($1*2.0/n):(mb_S($1*w,Ep/e)) w l lw 2 title "MB-6",\

set term wxt 1
set tics out nomirror
set title "Electron and Lattice energy distributions, 300K"
set ylabel "probablity"
set xlabel "e/E_f"
set boxwidth w/e
p "Electron_Velocity/5000.txt" u (bin(($5/e) -1.0,w/e)):(100.0/n) smooth freq w boxes title "500 fs",\
"" u ($1*w/e):(mb($1*w,Ee/e)) w l lw 2 title "MB",\
#"E_distrib" u (bin($2-e,w)):(100.0/n) smooth freq w boxes title "0 fs",\



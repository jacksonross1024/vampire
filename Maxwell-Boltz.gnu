#!/bin/bash

Ep = 300.0
Ee = 300.0
kB = 1.38e-3
e  = 105.59
n  = 8788.0
w  = 1.0/n

mb(x,t) = 0.15*sqrt(2.0/3.141592)*x*x*exp(-0.5*x*x/(t*t))/(t*t*t)
mb_e(x,t) = sqrt(2.0/3.141592)*x*x*exp(-1.0*x/t)/(t**1.5)
mb_TH(x,t) = (1.0/4.0)*sqrt(1*x/3.141592)*((1/t)**1.5)*exp(-x/t)
mb_S(x,t) = sqrt(3.141592)*(3.141592/16.0)*(x**2.0)*exp(-0.5*x/t)/(t**1.5)

bin(x,w) = w*floor(x/w)

#set term aqua 0
set term pngcairo font "arial, 14" size 1024,768
set output "MB-Distrib-CeTe-3ps-scaled.png"
set tics out nomirror
set title "Maxwell-Boltzmann Distribution C_e(T_e) Ni 3ps"
set ylabel "probablity"
set xlabel "e/E_f"
set boxwidth w
set xrange [1:1.1]
set auto y

p "0.5G/Atom_Energy/30000" u (bin($2/e,w)):(10.0/n) smooth freq w boxes title "0 ps",\
"" u (1.0+$1*w):(mb_S($1*w,300.0*kB/e)) w l lw 2 title "300K",\
"0.5G/Electron_Velocity/30000.txt" u (bin($2/e,w)):(10.0/n) smooth freq w boxes title "3 ps",\
"" u (1.0+$1*w):((1.0/sqrt(3.141592))*mb_S($1*w, 510.0*kB/e)) w l lw 2 title "510K",\



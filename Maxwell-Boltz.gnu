#!/bin/bash

T0 = 600
T1 = 1000
kB = 1.38e-3
e  = 105.59
n  = 8788.0
w  = 1.0/n
t0 = T0*kB/e
t1 = T1*kB/e

mb(x,t) = (0.01/16.0)*(x**2.0)*exp(-0.5*x/t)*(1.0/t)**3.0
mb_e(x,t) = sqrt(2.0/3.141592)*((x+3*t)**2)*exp(-1.0*(x+3*t)/t)/(t**1.5)
mb_TH(x,t) = (1.0/4.0)*sqrt(1*x/3.141592)*((1/t)**1.5)*exp(-x/t)
mb_S(x,t) = (0.01/16.0)*((x+(3*t))**2.0)*exp(-0.5*(x+(3*t))/t)*(1.0/t)**3.0

bin(x,w) = w*floor(x/w)

set term aqua 1
#set term pngcairo font "arial, 14" size 1024,768
#set output "MB-Electron-Lattice-Distrib100fs.png"
set tics out nomirror
set notitle #"Ideal Electron Gas, Electron and Lattice distributions, 200 fs"
set ylabel "probablity"
set xlabel "e/E_f"
set boxwidth w
set style fill solid 0.8 noborder
set xrange [1.0-3.0*t1:1.2]
set yrange [0:0.8]
#set auto x
#set auto y
set title "Classical Electron Gas, 20 mJ/cm^2"
set key top right notitle


p  "3G/Electron_Velocity/1000.csv" u (bin($2/e,w)):(100.0/n) smooth freq w boxes lc rgb "#20547C" title "100 fs",\
"3G/Electron_Velocity/6000.csv" u (bin($2/e,w)):(100.0/n) smooth freq w boxes lc rgb "#785F3D" title "600 fs",\
"" u (1.0-(3*t0)+$1*w):(mb_S(-3*t0+$1*w, t0)) w l lw 2 title sprintf("MB %.0fK", T0),\
"" u (1.0-(3*t1)+$1*w):(mb_S(-3*t1+$1*w, t1)) w l lw 2 title sprintf("MB %.0fK", T1),\

#"5G/Electron_Velocity/30000.txt" u (bin($2/e +1.0,w)-1.0):(100.0/n) smooth freq w boxes title "0.25e4",\
#"1G/Electron_Velocity/0.txt" u (bin($2/e +1.0,w)-1.0):(100.0/n) smooth freq w boxes title "0 ps",\

#"2G/Electron_Velocity/30000.txt" u (bin($2/e +1.0,w)-1.0):(100.0/n) smooth freq w boxes title "0.25e3",\

#"1G/Atom_Energy/30000" u (bin($2/e,w)):(1.0/n) smooth freq w boxes title "Lattice"

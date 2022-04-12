#!/bin/bash

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

set style line 1 pt 7 ps 1.2 lt 1 lc  rgb "#04233B" lw 2
set style line 2 pt 9 ps 1.6 lt 1 lc  rgb "#293E4D" lw 2
set style line 3 pt 5 ps 1.2 lt 1 lc  rgb "#20547C" lw 2
set style line 4 pt 11 ps 1.6 lt 1 lc rgb "#77B1DE" lw 2
set style line 5 pt 13 ps 1.6 lt 1 lc rgb "#EFBD79" lw 2
set style line 6 pt 15 ps 1.6 lt 1 lc rgb "#C07F28" lw 2
set style line 7 pt 9 ps 1.6 lt 1 lc  rgb "#785F3D" lw 2
set style line 8 pt 11 ps 1.6 lt 1 lc rgb "#5B3605" lw 2

set linetype 2 dt 2
#Set additional styles for dashed/dotted lines
set style line 100 pt 1 ps 1.2 lt 0 lc rgb "gray30" lw 2
set style line 101 pt 9 ps 1.4 lt 2 lc rgb "black" lw 2

#set term pngcairo font "arial, 14" size 1024,768 enhanced
#set output "TeTp-TTMe-TTLe-112G-3-5Tauee-hbaromega.png"
set term wxt 0
set tics in nomirror

set key top right
set style increment user


#===================================
#new plot for temperature/energy data
#===================================

set ylabel "x surface flux (counts)"
set xlabel "time (ps)"
set auto x
set auto y
set title "Free Electron Gas Model, 7.5V"
set key inside

p "1G/mean_data.csv" u ($1/1000):14 w l ls 1 title "1G",\
"2G/mean_data.csv" u ($1/1000):14 w l ls 2 title "2G",\
"3G/mean_data.csv" u ($1/1000):14 w l ls 3 title "3G",\
"5G/mean_data.csv" u ($1/1000):14 w l ls 4 title "5G"

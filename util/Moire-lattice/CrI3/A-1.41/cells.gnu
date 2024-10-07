

set palette defined (  0 'blue', 1 'white',  2 'red')

set term pngcairo 
set output "cells.png"

set ylabel "y position (nm)"
set xlabel "x position (nm)" 

set style fill solid noborder 



plot for[i=0:101*101:1]  "microcells.txt" every ::2 u (i%100):(floor(i/100)):(1.0):(1.0):(column(4+4*i)*column(5+4*i)) w boxxy lc palette notitle 


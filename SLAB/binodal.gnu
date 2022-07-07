#!/bin/gnuplot
set terminal pdfcairo enhanced color lw 3 size 4,3 font 'Arial-Bold'
set output 'Binodal.pdf'
set xlabel '{/*1.9{/Symbol:Italic:Bold r}}'
set ylabel '{/*1.9{/Arial-Italic:Bold T}}'

fg(x)=rhoc-A*(x-Tc)-B*(Tc-x)**0.32
fl(x)=rhoc-A*(x-Tc)+B*(Tc-x)**0.32
rhoc=0.22
A=0.21
Tc=1.17
B=-0.5

f(x,y) = y==0 ? fg(x) : fl(x)
fit f(x,y) "rho.txt" using 1:-2:2:(1) via rhoc,A,Tc,B

set xrange[0:0.9]
set yrange[0.65:1.20]

plot [0:0.9] [0.65:1.171] 'sq' u (rhoc-A*($1-Tc)+B*(Tc-$1)**0.32):($1) w l lw 0.5 lc rgb "black" notitle,  'sq' u (rhoc-A*($1-Tc)-B*(Tc-$1)**0.32):($1) w l lw 0.5 lc rgb "black" notitle,\
'rho.txt' u 2:1 pt 8 ps 1.0 lc rgb "blue" notitle,\
"<echo '0.313496 1.17009'" pt 8 ps 1.0 lc rgb "blue" notitle

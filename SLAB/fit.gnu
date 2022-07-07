#!/bin/gnuplot
set terminal pdfcairo enhanced color lw 3 size 4,3 font 'Arial-Bold'
set output 'HistogramFit.pdf'
set xlabel '{/*1.9{/Arial-Italic:Bold z}}'
set ylabel '{/*1.9{/Symbol:Italic:Bold r}}'
set yrange[0:0.9]
set xrange[0:0]
set multiplot


set xtic font ",18"
set ytic font ",18"
set key font ",18"
set key font ",18"

f(x) = a + c*tanh(d*(x-e))

a=0.4128
c=-0.4129
d=1.378
e=3.05
fit [0:75] [*:*]  f(x) 'T0.65/histograms.txt' using 1:2 via a,c,d,e

set print "output.txt"
print 0.65,a+c,a-c
set print

fit [0:75] [*:*]  f(x) 'T0.70/histograms.txt' using 1:2 via a,c,d,e

set print "output.txt" append
print 0.70,a+c,a-c
set print

fit [0:75] [*:*]  f(x) 'T0.75/histograms.txt' using 1:2 via a,c,d,e

set print "output.txt" append
print 0.75,a+c,a-c
set print

fit [0:75] [*:*]  f(x) 'T0.80/histograms.txt' using 1:2 via a,c,d,e

set print "output.txt" append
print 0.80,a+c,a-c
set print

fit [0:75] [*:*]  f(x) 'T0.85/histograms.txt' using 1:2 via a,c,d,e

set print "output.txt" append
print 0.85,a+c,a-c
set print

fit [0:75] [*:*]  f(x) 'T0.90/histograms.txt' using 1:2 via a,c,d,e

set print "output.txt" append
print 0.90,a+c,a-c
set print

fit [0:75] [*:*]  f(x) 'T0.95/histograms.txt' using 1:2 via a,c,d,e

set print "output.txt" append
print 0.95,a+c,a-c
set print

fit [0:75] [*:*]  f(x) 'T1.00/histograms.txt' using 1:2 via a,c,d,e

set print "output.txt" append
print 1.0,a+c,a-c
set print

fit [0:75] [*:*]  f(x) 'T1.05/histograms.txt' using 1:2 via a,c,d,e

set print "output.txt" append
print 1.05,a+c,a-c
set print

fit [0:75] [*:*]  f(x) 'T1.10/histograms.txt' using 1:2 via a,c,d,e

set print "output.txt" append
print 1.10,a+c,a-c
set print

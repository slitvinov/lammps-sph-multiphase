set term push
set term postscript eps enhanced color
set output sprintf("%s.eps", "$0")
replot

set terminal png
set output sprintf("%s.png", "$0")
replot
set output

set term pop
replot
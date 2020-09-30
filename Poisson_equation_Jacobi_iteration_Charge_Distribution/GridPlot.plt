set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "GridPlot.pdf"
set autoscale
set grid
set label
set title "Grid Evolution"
set xlabel "X" 
set ylabel "Y" 
set zlabel "Phi(x,y)" 
set key autotitle
splot "GridStep00.dat" using 1:2:3 with points title "Starting grid" lc "red"
splot "GridStepND.dat" using 1:2:3 with points title "Final grid" lc "blue"

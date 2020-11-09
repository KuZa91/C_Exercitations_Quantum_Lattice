set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "AvgCorrelator.pdf"
set autoscale
set grid
set label
set title "Average Correlator among Samples"
set xlabel "Tau" 
set ylabel "<Phi(Tau)Phi(0)>" 
set key autotitle
plot "avg_correlator.dat" using 1:2 with points title "Correlator Values" lc "red"
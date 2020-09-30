set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "IsingMagnetization.pdf"
set autoscale
set title "Average Spin vs Iteration"
unset key
set xlabel "Iteration [n]"
set ylabel "<S>"
plot "State_Variables.dat" using 1:2 with lines title "<S>(n)"

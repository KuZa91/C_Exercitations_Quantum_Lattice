set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "WavePlot.pdf"
set autoscale
set grid
set label
set title "Wave Evolution"
set xlabel "X" 
set ylabel "Real(Waveform(x))" 
set key autotitle
plot "WaveEvol00.dat" using 1:($2**2 + $3**2) with points title "Starting Wave" lc "red"
plot "WaveEvolND.dat" using 1:($2**2 + $3**2) with points title "Final Wave" lc "blue"

reset

# png
set terminal pngcairo size 1080,1920 enhanced font 'Times New Roman,8'

# color definitions
set border linewidth 2
set style line 1 lc rgb 'grey30' ps 0 lt 1 lw 2
set style line 2 lc rgb 'grey70' lt 1 lw 2
set style fill solid 1.0 border rgb 'grey30'

set xrange [500:2000]
set xtics nomirror 
set ytics nomirror
set datafile separator ','

set output 'Output/Outputs.png'

### Start multiplot (2x2 layout)
set multiplot layout 4,1 rowsfirst
# --- GRAPH a
#set label 1 'Energy' at graph 0.92,0.9 font ',8'
set ylabel 'Energy'
set xlabel 'Step'
set title "Energy vs Time"
plot 'Data.csv' using 1:6 with lines title 'Box A', '' using 1:7 with lines title 'Box B'

# --- GRAPH b
#set label 1 'b' at graph 0.92,0.9 font ',8'
set ylabel 'Density'
set xlabel 'Step'
set title "Density vs Time"
plot 'Data.csv' using 1:8 with lines title 'Box A', '' using 1:9 with lines title 'Box B'

# --- GRAPH c
#set label 1 'c' at graph 0.92,0.9 font ',8'
set ylabel 'Population of Particles'
set xlabel 'Step'
set yrange [100:150]
set title "Population vs Time"
plot 'Data.csv' using 1:4 with lines title 'Box A', '' using 1:5 with lines title 'Box B'
unset yrange

# --- GRAPH d
#set label 1 'd' at graph 0.92,0.9 font ',8'
set ylabel 'Volume'
set xlabel 'Step'
set title "Volume vs Time"
plot 'Data.csv' using 1:($2**3) with lines title 'Box A', '' using 1:($3**3) with lines title 'Box B'
### End multiplot
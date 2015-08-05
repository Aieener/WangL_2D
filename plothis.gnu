# plothis.gnu
set terminal png
set term png truecolor
set output "Histo.png"
set title "Histogram"
set xlabel "Ni"
set ylabel "Freq(Ni)"
binwidth=10
bin(x,width)=width*floor(x/width)
set style fill transparent solid 0.3333 border -1
# Add light grid to plot
set style line 102 lc rgb '#d6d7d9' lt 0 lw 1
set grid back ls 102

# xyborder.cfg
set style line 101 lc rgb '#808080' lt 1 lw 1
set border 3 front ls 101
set tics nomirror out scale 0.75
set format '%g'

plot 'datahis.dat' using (bin($5,binwidth)):(1) smooth freq with boxes title "Nv",'dataplot.dat' using (bin($4,binwidth)):(1) smooth freq with boxes title "Nh"

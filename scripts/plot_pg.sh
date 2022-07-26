base=${1:-"EFGHJP"}
run=${2:-"run100"}
epoch=${3:-"01"}

base2=`echo $base | sed 's/_/-/g'`
echo $base2
gnuplot<<EOF
set term pdf
set out "figs/${base}_${run}_${epoch}.pdf"

set xlabel "Genotype"
set ylabel "Phenotype"
set grid
set key title 'generation'
set key bottom

set xrange [-0.09:1.09]
set yrange [-0.09:1.09]
set xtics 0.1
set ytics 0.1

set pointsize 0.5

set title "Model: ${base2} ${epoch}"
p "< awk '/Data/' proj/${base}_${run}_${epoch}_001.dat" u 4:5 t '1'  w p pt 7 lt rgb "#0000FF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_010.dat" u 4:5 t '10' w p pt 7 lt rgb "#6020FF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_025.dat" u 4:5 t '25' w p pt 7 lt rgb "#8020FF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_050.dat" u 4:5 t '50' w p pt 7 lt rgb "#A020FF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_100.dat" u 4:5 t '100' w p pt 7 lt rgb "#B020FF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_150.dat" u 4:5 t '150' w p pt 7 lt rgb "#C020FF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_200.dat" u 4:5 t '200' w p pt 7 lt rgb "#FF20FF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_001.dat" u 2:3 t '' w p pt 6 lt rgb "#0000FF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_010.dat" u 2:3 t '' w p pt 6 lt rgb "#602AFF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_025.dat" u 2:3 t '' w p pt 6 lt rgb "#8020FF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_050.dat" u 2:3 t '' w p pt 6 lt rgb "#A020FF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_100.dat" u 2:3 t '' w p pt 6 lt rgb "#B020FF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_150.dat" u 2:3 t '' w p pt 6 lt rgb "#C020FF" \
, "< awk '/Data/' proj/${base}_${run}_${epoch}_200.dat" u 2:3 t '' w p pt 6 lt rgb "#FF20FF" 

set out "figs/${base}_${run}_${epoch}_nov.pdf"
set xlabel "Genotype"
set ylabel "Phenotype"
set grid
set key bottom

set pointsize 0.2
set key title "Novel Environment"
set title "Model: ${base2} ${epoch}"
p "<awk '/AVESD1/' proj/${base}_${run}_${epoch}_*.dat" u 3:4:5:6 t '' w xyerrorlines lt rgb "#1010FF" 

reset
set out "figs/${base}_${run}_${epoch}_anc.pdf"
set xlabel "Genotype"
set ylabel "Phenotype"
set grid
set key bottom
set key title "Ancestral Environment"
set xrange [-0.09:1.09]
set yrange [-0.09:1.09]
set xtics 0.1
set ytics 0.1

set pointsize 0.5
set title "Model: ${base2} ${epoch}"
p "<awk '/AVESD0/' proj/${base}_${run}_${epoch}_*.dat" u 3:4:5:6 t '' w xyerrorlines lt rgb "#FF1010" 

reset
set out "figs/${base}_${run}_${epoch}_both.pdf"
set xlabel "Genotype"
set ylabel "Phenotype"
set grid
set key bottom
set xrange [-0.09:1.09]
set yrange [-0.09:1.09]
set xtics 0.1
set ytics 0.1

set pointsize 0.5
set title "Model: ${base2} ${epoch}"
p "<awk '/AVESD1/' proj/${base}_${run}_${epoch}_*.dat" u 3:4:5:6 t 'Novel Environment' w xyerrorlines lt rgb "#1010FF" \
, "<awk '/AVESD0/' proj/${base}_${run}_${epoch}_*.dat" u 3:4:5:6 t 'Ancestral Environment' w xyerrorlines lt rgb "#FF1010" 
EOF


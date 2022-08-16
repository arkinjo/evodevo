
base=${1:-"EFGHJP"}
run=${2:-"run100"}
epoch=${3:-"01"}

zsh plot_pg.sh $base "${run}" $epoch 

base2=`echo $base | sed 's/_/-/g'`
echo $base2
filename=${base}_tmpanim.gnu
cat > ${filename}<<EOF
set term pdf
set out "figs/${base}_${run}_${epoch}_anim.pdf"
EOF
for i in {001..200}; do
cat<<EOF
reset
set xlabel "Genotype"
set ylabel "Phenotype"
set grid
set key bottom

set xrange [-0.09:1.09]
set yrange [-0.09:1.09]
set xtics 0.1
set ytics 0.1
#set pointsize 0.5

set title "Model: ${base2}; generation ${i}"
p "< awk '/Data/' proj/${base}_${run}_${epoch}_${i}.dat" u 4:5 t 'Novel Environment'  w p pt 7 lt rgb "#1010FF" \
, '' u 2:3 t 'Ancestral Environment'  w p pt 6 lt rgb "#FF1010" \
, "< awk '/AVESD1/' proj/${base}_${run}_${epoch}_${i}.dat" u 3:4:5:6 t ''  w xyerrorbars pt 1 lt rgb "#FF10FF" \
, "< awk '/AVESD0/' proj/${base}_${run}_${epoch}_${i}.dat" u 3:4:5:6 t ''  w xyerrorbars pt 1 lt rgb "#10FFFF"
EOF
done >> ${filename}
gnuplot ${filename}
rm -f ${filename}

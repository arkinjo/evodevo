base=${1:-"EFGHJP"}
#base=EFGHJP
base2=`echo $base | sed 's/_/-/g'`
gnuplot<<EOF
set term pdf enhanced
set encoding iso_8859_1

set title "model: ${base2}"
set ylabel "{/Symbol Dd}p (projection on {/Symbol d}{/Times-Bold \352})"
set key title "{/Symbol d}{/Times-Bold e} (%)"
set key bottom

set grid
set pointsize 0.5
set xrange [-10:10]
set yrange [-10:10]
set out "covprjpe_${base}_5_50.pdf"
set xlabel "{/Symbol Dd}e (projection on {/Symbol d}{/Times-Bold \352}^T<{/Symbol Dd}p{/Symbol Dd}e^T>)"
stats "< awk /Prj/ per/${base}_10_0.05_*.dat " u 4:6 name "e05" 
stats "< awk /Prj/ per/${base}_100_0.05_*.dat " u 4:6 name "e50"

p "< awk /Prj/ per/${base}_10_0.05_*.dat " u 4:6 t "5" lt 1 \
, "< awk /Prj/ per/${base}_100_0.05_*.dat " u 4:6 t "50" lt 2 \
, e05_slope*x t sprintf("%4.2f x (r=%4.2f)", e05_slope, e05_correlation) lt 1 \
, e50_slope*x t sprintf("%4.2f x (r=%4.2f)", e50_slope, e50_correlation) lt 2 

set xrange [-2:2]
set out "covprjpG_${base}_5_50.pdf"
set xlabel "{/Symbol Dd}{/Times-Bold G} (projection on {/Symbol d}{/Times-Bold \352}^T<{/Symbol Dd}p{/Symbol Dd}{/Times-Bold G}^T>)"
stats "< awk /Prj/ per/${base}_10_0.05_*.dat " u 5:6 name "g05"
stats "< awk /Prj/ per/${base}_100_0.05_*.dat " u 5:6 name "g50"

p "< awk /Prj/ per/${base}_10_0.05_*.dat " u 5:6 t "5" \
, "< awk /Prj/ per/${base}_100_0.05_*.dat " u 5:6 t "50"\
, g05_slope*x t sprintf("%4.2f x (r=%4.2f)", g05_slope, g05_correlation) lt 1\
, g50_slope*x t sprintf("%4.2f x (r=%4.2f)", g50_slope, g50_correlation) lt 2
EOF

###############################################################################

for denv in 2 {10..100..10}; do
    xx=$(($denv / 2))
cat<<EOF
set term pdf enhanced
set encoding iso_8859_1

set title "model: ${base2}"
set ylabel "{/Symbol Dd}p (projection on {/Symbol d}{/Times-Bold \352})"
set key title "{/Symbol d}{/Times-Bold e} = ${xx}%"
set key bottom

set grid
set pointsize 0.5
set xrange [-10:10]
set yrange [-10:10]
set xlabel "{/Symbol Dd}e (projection on {/Symbol d}{/Times-Bold \352}^T<{/Symbol Dd}p{/Symbol Dd}e^T>)"
set out "figs/covprjpe_${base}_${denv}.pdf"
stats "< awk /Prj/ per/${base}_${denv}_0.05_*.dat " u 4:6 name "es"
p "< awk /Prj/ per/${base}_${denv}_0.05_*.dat " u 4:6 t "" lt 1 \
, es_slope*x t sprintf("%4.2f x (r=%4.2f)", es_slope, es_correlation) lt 1

###########################
set xrange [-2:2]
set xlabel "{/Symbol Dd}{/Times-Bold G} (projection on {/Symbol d}{/Times-Bold \352}^T<{/Symbol Dd}p{/Symbol Dd}{/Times-Bold G}^T>)"
set out "figs/covprjpG_${base}_${denv}.pdf"
stats "< awk /Prj/ per/${base}_${denv}_0.05_*.dat " u 5:6 name "gs"
p "< awk /Prj/ per/${base}_${denv}_0.05_*.dat " u 5:6 t "" \
, gs_slope*x t sprintf("%4.2f x (r=%4.2f)", gs_slope, gs_correlation) lt 1

EOF
done | gnuplot

pdfjoin -o covprjpe_${base}_anim.pdf figs/covprjpe_${base}_2.pdf  figs/covprjpe_${base}_{10..100..10}.pdf
pdfjoin -o covprjpG_${base}_anim.pdf figs/covprjpG_${base}_2.pdf  figs/covprjpG_${base}_{10..100..10}.pdf


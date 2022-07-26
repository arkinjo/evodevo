for base in EFGHJP _FGHJP E_G__P EFGHJ_; do
    base2=`echo $base | sed 's/_/-/g'`
    for i in {01..20}; do
	gnuplot<<EOF
set term pdf
set xlabel "Singular component"
set ylabel "Cumulative singular values (%)"
set key title "{/Symbol d}{/Times-Bold e} (%)"
set key bottom
set pointsize 0.5
set yrange [0:100]

set out "figs/${base}_svcum${i}.pdf"
set title "${base2}"
p "< grep SVal per/${base}_2_0.05_${i}.dat" u (\$2+1):(\$5*100) t '1' w lp \
, "< grep SVal per/${base}_10_0.05_${i}.dat" u (\$2+1):(\$5*100) t '5' w lp \
, "< grep SVal per/${base}_20_0.05_${i}.dat" u (\$2+1):(\$5*100) t '10' w lp \
, "< grep SVal per/${base}_30_0.05_${i}.dat" u (\$2+1):(\$5*100) t '15' w lp \
, "< grep SVal per/${base}_40_0.05_${i}.dat" u (\$2+1):(\$5*100) t '20' w lp \
, "< grep SVal per/${base}_50_0.05_${i}.dat" u (\$2+1):(\$5*100) t '25' w lp \
, "< grep SVal per/${base}_60_0.05_${i}.dat" u (\$2+1):(\$5*100) t '30' w lp \
, "< grep SVal per/${base}_70_0.05_${i}.dat" u (\$2+1):(\$5*100) t '35' w lp \
, "< grep SVal per/${base}_80_0.05_${i}.dat" u (\$2+1):(\$5*100) t '40' w lp \
, "< grep SVal per/${base}_90_0.05_${i}.dat" u (\$2+1):(\$5*100) t '45' w lp \
, "< grep SVal per/${base}_100_0.05_${i}.dat" u (\$2+1):(\$5*100) t '50' w lp

EOF

    done
    pdfjoin -o figs/${base}_svcum.pdf figs/${base}_svcum{01..20}.pdf
done
    

base=${1:-"EFGHJP"}
ccamode=${2:-"2"} #Mode of cross covariance analysis: 0: ancestral environment, 1: novel environment, 2: difference between environments


RANDOM=1339
NIND=20
NOISE=0.05

## cross-covariance analysis
pegcov=../GitHub/evodevo/crosscov/crosscov
#


# perturb and compute cross-covariance for the 1st generations.
for denv in 2 {10..100..10}; do
    for i in {01..${NIND}}; do
	zsh perturb.sh $base ${denv} ${NOISE} $i $RANDOM
	$pegcov -jsonin pops/${base}_perturb${denv}_${NOISE}_01_001.json.gz \
		-mode=$ccamode -right=e \
		> per/${base}_${denv}_${NOISE}_${ccamode}_${i}_pe.dat
	$pegcov -jsonin pops/${base}_perturb${denv}_${NOISE}_01_001.json.gz \
		-mode=$ccamode -right=G \
		> per/${base}_${denv}_${NOISE}_${ccamode}_${i}_pG.dat

    done
done

mkdir LSVec1 LProj1

# extracting data
for cc in pe pG; do
    for denv in 2 {10..100..10}; do
	for i in {01..${NIND}}; do
	    echo -n "${denv}"
	    grep FN2 per/${base}_${denv}_${NOISE}_${ccamode}_${i}_${cc}.dat | awk '{printf "\t%e", $2}'
	    echo ''
	done
    done > ${base}_${ccamode}_${cc}_denv22.dat

    # extracting data: Alignment between env. change and PCaxis 
    for denv in 2 {10..100..10}; do
	for i in {01..${NIND}}; do
	    echo -n $denv
	    grep "Ali" per/${base}_${denv}_${NOISE}_${ccamode}_${i}_${cc}.dat | awk '$2==0 {printf "\t%f\n", $3}'
	done
    done > ${base}_${ccamode}_${cc}_ali.dat

    # # This "Sigma" is a bit complicated. ||environmental change x cross-cov||
    # for denv in 2 {10..100..10}; do
    # 	for i in {01..${NIND}}; do
    # 	    echo -n $denv
    # 	    grep "Sigma" per/${base}_${denv}_${NOISE}_${i}_${cc}.dat | awk '{printf "\t%e\n", $2}'
    # 	done
    # done > ${base}_${cc}_sigma.dat

    # for denv in 2 {10..100..10}; do
    # 	for i in {01..${NIND}}; do
    # 	    echo -n $denv
    # 	    fgrep Pdenv per/${base}_${denv}_${NOISE}_${i}_${cc}.dat | awk '{printf "\t%e\t%e\t%e\n", $2, $3, $4}'
    # 	done
    # done > ${base}_${cc}_pheno.dat

    # 1st singular values
    for denv in 2 {10..100..10}; do
	for i in {01..${NIND}}; do
	    echo -n $denv
	    fgrep SVal per/${base}_${denv}_${NOISE}_${ccamode}_${i}_${cc}.dat | awk '$2==0 {printf "\t%e\t%e\t%e\n", $3, $4, $5}'
	done
    done > ${base}_${ccamode}_${cc}_sval1.dat

    # 2nd singular values
    for denv in 2 {10..100..10}; do
	for i in {01..${NIND}}; do
	    echo -n $denv
	    fgrep SVal per/${base}_${denv}_${NOISE}_${ccamode}_${i}_${cc}.dat | awk '$2==1 {printf "\t%e\t%e\t%e\n", $3, $4, $5}'
	done
    done > ${base}_${ccamode}_${cc}_sval2.dat

    # cumulative singular values except 1st SV. ("everything else")
    for denv in 2 {10..100..10}; do
	for i in {01..${NIND}}; do
	    echo -n $denv
	    fgrep SVal per/${base}_${denv}_${NOISE}_${ccamode}_${i}_${cc}.dat | awk '$2==0 {printf "\t%e\n", $3*(1.0/$4 - 1.0)}'
	done
    done > ${base}_${ccamode}_${cc}_svals.dat

done

# First singular vector elements
for denv in 2 {10..100..10}; do
    for i in {01..${NIND}}; do
	fgrep LSVec1 per/${base}_${denv}_${NOISE}_${ccamode}_${i}_pe.dat > tmp0.dat
	fgrep LSVec1 per/${base}_${denv}_${NOISE}_${ccamode}_${i}_pG.dat > tmp1.dat
	paste tmp0.dat tmp1.dat | awk '{printf "%d\t%e\t%e\n", $2, $3, $6}' > LSVec1/${base}_${denv}_${NOISE}_${ccamode}_${i}.lsvec1
    done 
done
rm -f tmp0.dat tmp1.dat

# Phenotype projected to singular vector
for denv in 2 {10..100..10}; do
    for i in {01..${NIND}}; do
    fgrep LProj1 per/${base}_${denv}_${NOISE}_${ccamode}_${i}_pe.dat > tmp0.dat
    fgrep LProj1 per/${base}_${denv}_${NOISE}_${ccamode}_${i}_pG.dat > tmp1.dat
    paste tmp0.dat tmp1.dat | awk '{printf "%d\t%e\t%e\n", $2, $3, $6}' > LProj1/${base}_${denv}_${NOISE}_${ccamode}_${i}.lproj1
    done
done
rm -f tmp0.dat tmp1.dat
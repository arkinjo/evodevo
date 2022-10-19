base=${1:-"EFGHJP"}

RANDOM=1339
NIND=20
NOISE=0.05

## cross-covariance analysis
pegcov=~/Documents/GitHub/evodevo/pegcov/pegcov
#

# perturb and compute cross-covariance for the 1st generations.
for denv in 2 {10..100..10}; do
    for i in {01..${NIND}}; do
	zsh perturb.sh $base ${denv} ${NOISE} $i $RANDOM
	$pegcov -jsonin pops/${base}_perturb${denv}_${NOISE}_01_001.json.gz \
		-p=2 -eg=2 \
		> per/${base}_${denv}_${NOISE}_${i}.dat
    done
done

# extracting data
for denv in 2 {10..100..10}; do
    for i in {01..${NIND}}; do
	echo -n "${denv}"
	grep Ddp per/${base}_${denv}_${NOISE}_${i}.dat | awk '{printf "\t%e", $2}'
	echo ''
    done
done > ${base}_denv22.dat

# extracting data: Alignment between env. change and PCaxis 
for denv in 2 {10..100..10}; do
    for i in {01..${NIND}}; do
	echo -n $denv
	grep "Ali" per/${base}_${denv}_${NOISE}_${i}.dat | awk '$2==0 {printf "\t%f\n", $3}'
    done
done > ${base}_ali.dat

# This "Sigma" is a bit complicated. ||environmental change x cross-cov||
for denv in 2 {10..100..10}; do
    for i in {01..${NIND}}; do
	echo -n $denv
	grep "Sigma" per/${base}_${denv}_${NOISE}_${i}.dat | awk '{printf "\t%e\n", $2}'
    done
done > ${base}_sigma.dat

for denv in 2 {10..100..10}; do
    for i in {01..${NIND}}; do
	echo -n $denv
	fgrep Pdenv per/${base}_${denv}_${NOISE}_${i}.dat | awk '{printf "\t%e\t%e\t%e\n", $2, $3, $4}' #Cross covariance between phenotype and environment
    done
done > ${base}_pheno.dat

# 1st singular values
for denv in 2 {10..100..10}; do
    for i in {01..${NIND}}; do
	echo -n $denv
	fgrep SVal per/${base}_${denv}_${NOISE}_${i}.dat | awk '$2==0 {printf "\t%e\t%e\t%e\n", $3, $4, $5}' #0+1th singular value, prop, cum prop
    done
done > ${base}_sval1.dat

# 2nd singular values
for denv in 2 {10..100..10}; do
    for i in {01..${NIND}}; do
	echo -n $denv
	fgrep SVal per/${base}_${denv}_${NOISE}_${i}.dat | awk '$2==1 {printf "\t%e\t%e\t%e\n", $3, $4, $5}'
    done
done > ${base}_sval2.dat

# cumulative singular values except 1st SV. ("everything else")
for denv in 2 {10..100..10}; do
    for i in {01..${NIND}}; do
	echo -n $denv
	fgrep SVal per/${base}_${denv}_${NOISE}_${i}.dat | awk '$2==0 {printf "\t%e\n", $3*(1.0/$4 - 1.0)}'
    done
done > ${base}_svals.dat

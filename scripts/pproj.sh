base=${1:-"EFGHJP"}
#ccamode=${2:-"2"} #Mode of cross covariance analysis: 0: ancestral environment, 1: novel environment, 2: difference between environments


RANDOM=20101019
NIND=20
NOISE=0.05

## project phenotypes
pproj=../GitHub/evodevo/pproj/pproj
#


# perturb and compute cross-covariance for the 1st generations.
for denv in 2 {10..100..10}; do
    for i in {01..${NIND}}; do
	zsh perturb.sh $base ${denv} ${NOISE} $i $RANDOM
	$pproj -jsongzin pops/${base}_perturb${denv}_${NOISE}_01_001.json.gz \
		-PG_File pproj/${base}_run${denv}_${epo} \
    done
done

#!/bin/zsh

export GOMAXPROCS=16

base=${1:-"EFGHJP"} # E_G__P 
DENV2=${2:-"100"}
RANDOM=${3:-"13531"}

EVODEVODIR=~/work/GitHub/evodevo/
train=${EVODEVODIR}/train/train
#ccphenv=${EVODEVODIR}/ccphenv/ccphenv
pgproj=${EVODEVODIR}/pgproj/pgproj


MAXPOP=200
NEPOCH1=40
NEPOCH2=10
NCELLS=1
DENV1=100

NOISE=0.05
MUT=0.001
NGEN=200
SEED1=${RANDOM}
SEEDCUE1=${RANDOM}
SEED2=${RANDOM}
SEEDCUE2=${RANDOM}
REF1=001
REF2=200
TAUF=1

# default values
if [ "${base}" = "E_G__P" ]; then
#    Single-layer with fat hidden layer.
    NGENES=600
    DENSITY_E=$((0.02/3))
    DENSITY_G=$((4*0.02/9))
    DENSITY_P=$((0.02/3))
else
    NGENES=200
    DENSITY_E=0.02
    DENSITY_G=0.02
    DENSITY_P=0.02
fi


if [ ${base[1]} = "_" ]; then e=F; else e=T; fi
if [ ${base[2]} = "_" ]; then f=F; else f=T; fi
if [ ${base[4]} = "_" ]; then h=F; else h=T; fi
if [ ${base[5]} = "_" ]; then j=F; else j=T; fi
if [ ${base[6]} = "_" ]; then p=F; else p=T; fi

echo $base ${e}${f}G${h}${j}${p}

if [ ${DENV2} -eq 100 ]; then
    $train -maxpop=${MAXPOP} -ncells=${NCELLS} -nepoch=${NEPOCH1} \
	   -traj_file=traj/${base}_train.traj \
	   -jsonout=json/${base}_train.json \
	   -denv=${DENV1} -noise=${NOISE} -mut=${MUT} \
	   -cue=${e} -layerF=${f} -layerH=${h} -layerJ=${j} \
	   -pfback=${p} \
	   -seed=${SEED1} -seed_cue=${SEEDCUE1} \
	   -ngenes=${NGENES} -dE=${DENSITY_E} -dG=${DENSITY_G} -dP=${DENSITY_P} \
	   -tauF=${TAUF} \
	   > /dev/null
fi

$train -test=true -nepoch=${NEPOCH2} -maxpop=${MAXPOP} -ncells=${NCELLS} \
       -traj_file=traj/${base}_run${DENV2}.traj \
       -jsonin=json/${base}_train.json \
       -jsonout=pops/${base}_run${DENV2} \
       -denv=${DENV2} -noise=${NOISE} -mut=${MUT} \
       -seed=${SEED2} -seed_cue=${SEEDCUE2} \
       > /dev/null

for epo in {01..${NEPOCH2}}; do
    $pgproj -maxpop=${MAXPOP} -ncells=${NCELLS} -ngen=${NGEN} \
	    -ref1=pops/${base}_run${DENV2}_${epo}_${REF1}.json \
	    -ref2=pops/${base}_run${DENV2}_${epo}_${REF2}.json \
	    -jsonin=pops/${base}_run${DENV2}_${epo} \
	    -PG_file=proj/${base}_run${DENV2}_${epo} -env=false
done


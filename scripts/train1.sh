#!/bin/zsh

export GOMAXPROCS=16

base=${1:-"EFGHJP"} # E_G__P 
DENV2=${2:-"100"}
RANDOM=${3:-"13531"}

EVODEVODIR=~/Documents/GitHub/evodevo
train=${EVODEVODIR}/train/train
#ccphenv=${EVODEVODIR}/ccphenv/ccphenv
pgproj=${EVODEVODIR}/pgproj/pgproj


MAXPOP=200
MAXDEVSTEP=200
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
#    Single-layer with fat hidden layer. [NoHier]
    NGENES=600
    DENSITY_E=$((0.02/3))
    DENSITY_G=$((4*0.02/9))
    DENSITY_H=0.00
    DENSITY_P=$((0.02/3))
elif [ "${base}" = "EFGH__" ]; then
#   "Feedforward" deep model [NoDev]
    NGENES=200
    DENSITY_E=0.02
    DENSITY_G=0.02
    DENSITY_H=0.04
    DENSITY_P=0.02
    MAXDEVSTEP=1
elif [ "${base}" = "__G___" ]; then
#    Single fat layer "feedforward" model. [Null]   
    NGENES = 800
    DENSITY_E=0.00
    DENSITY_G=$((5*0.02/16))
    DENSITY_H=0.00
    DENSITY_P=$((0.02/4))
else
#   Everything else
    NGENES=200
    DENSITY_E=0.02
    DENSITY_G=0.02
    DENSITY_H=0.02
    DENSITY_P=0.02
fi


if [ ${base[1]} = "_" ]; then e=F; else e=T; fi
if [ ${base[2]} = "_" ]; then f=F; else f=T; fi
if [ ${base[4]} = "_" ]; then h=F; else h=T; fi
if [ ${base[5]} = "_" ]; then j=F; else j=T; fi
if [ ${base[6]} = "_" ]; then p=F; else p=T; fi

echo $base ${e}${f}G${h}${j}${p}

if [ ${DENV2} -eq 100 ]; then
    $train -maxpop=${MAXPOP} -maxdevstep=${MAXDEVSTEP} -ncells=${NCELLS} -nepoch=${NEPOCH1} \
	   -traj_file=traj/${base}_train.traj \
	   -jsongzout=json/${base}_train.json.gz \
	   -denv=${DENV1} -noise=${NOISE} -mut=${MUT} \
	   -cue=${e} -layerF=${f} -layerH=${h} -layerJ=${j} \
	   -pfback=${p} \
	   -seed=${SEED1} -seed_cue=${SEEDCUE1} \
	   -ngenes=${NGENES} -dE=${DENSITY_E} -dG=${DENSITY_G} -dH=${DENSITY_H} -dP=${DENSITY_P} \
	   -tauF=${TAUF} \
	   > /dev/null
fi

for denv in 10; do 
    $train -test=true -nepoch=${NEPOCH2} -maxpop=${MAXPOP} -ncells=${NCELLS} \
       -traj_file=traj/${base}_run${denv}.traj \
       -jsongzin=json/${base}_train.json.gz \
       -jsongzout=pops/${base}_run${denv} \
       -denv=${denv} -noise=${NOISE} -mut=${MUT} \
       -seed=${SEED2} -seed_cue=${SEEDCUE2} \
       > /dev/null

    for epo in {01..${NEPOCH2}}; do
        $pgproj -maxpop=${MAXPOP} -ncells=${NCELLS} -ngen=${NGEN} \
	        -ref1=pops/${base}_run${denv}_${epo}_${REF1}.json.gz \
	        -ref2=pops/${base}_run${denv}_${epo}_${REF2}.json.gz \
	        -jsongzin=pops/${base}_run${denv}_${epo} \
	        -PG_file=proj/${base}_run${denv}_${epo} -env=false
    done
done


#!/bin/zsh

export GOMAXPROCS=8
base=${1:-"EFGHJP"}
DENV=${2:-"0"}
NOISE=${3:-"0.00"}
IND=${4:-"01"}
RANDOM=${5:-"13531"}

EVODEVODIR=~/Documents/GitHub/evodevo
train=${EVODEVODIR}/train/train
#ccphenv=${EVODEVODIR}/ccphenv/ccphenv
pgproj=${EVODEVODIR}/pgproj/pgproj

export GOMAXPROCS=8
MAXPOP=200
NEPOCH1=40
NEPOCH2=1
NCELLS=1
MUT=0.001
#NOISE=0.1
NGEN=1
SEED1=${RANDOM}
SEEDCUE1=${RANDOM}
SEED2=${RANDOM}
SEEDCUE2=${RANDOM}
REF1=001
REF2=200
# run E F H J P
echo $base

$train -test=true -nepoch=${NEPOCH2} -maxpop=${MAXPOP} -ncells=${NCELLS} \
       -traj_file=traj/${base}_perturb${DENV}_${NOISE}_${IND}.traj \
       -jsongzin=json/${base}_train.json.gz \
       -jsongzout=pops/${base}_perturb${DENV}_${NOISE} \
       -seed=${SEED2} -seed_cue=${SEEDCUE2} \
       -ngen=${NGEN} -denv=${DENV} -noise=${NOISE} -mut=${MUT}\
       > gegege

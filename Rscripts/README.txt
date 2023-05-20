PGtraj.R generates phenotype-genotype plots and should be run in the directory containing the outputs of pgproj.go.

apr.R generates plots pertaining to plastic response in the first generation and uncovering of cryptic mutations. It should be run in the directory containing the outputs of pgproj.go.

ModelComp.R generates plots from SVD analysis. It should be run in the directory containing the outputs of check_perturb2.sh

lsvec1.R plots alignment between Pheno-Cue and Pheno-Geno axis. It should be run in the directory containing .lsvec1 files.

gvar.R generates plots pertaining to genetic variance and genetic bottleneck. It should be run in the directory containing the .gvar files.

traj.R generates other miscellaneous plots, in particular, those pertaining to the mismatch of the phenotype.

gvarANDpproj.R should only be run after running gvar.R and PGtraj.R. It generates plots of projected phenotype and genotype at the genetic bottleneck.

gvarANDerr.R should only be run after running gvar.R and traj.R. It generates plots of mismatch before and after the genetic bottleneck

*-train.R functions the same as it's above variants, but compares Full and NoTrain models instead of models of different features.
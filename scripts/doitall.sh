DENV2=100

mkdir pops traj proj json per cov figs

for base in EFGHJP _FGHJP E_G__P EFGHJ_; do
    (
	# train, test & project
	zsh train1.sh $base $DENV2

	# requires gnuplot
	#zsh make_animation.sh $base run${DENV2}
	#for i in {01..10}; do
    #	    zsh plot_pg.sh $base run${DENV2} $i
	#done
    #)
    # perturbation & check 1st generations
    for mode in {0..2} ; do 
        zsh check_perturb2.sh $base $mode #Loop over all modes
    done

    # requires gnuplot
    #zsh plot_covprj.sh $base 
	)
done

wait

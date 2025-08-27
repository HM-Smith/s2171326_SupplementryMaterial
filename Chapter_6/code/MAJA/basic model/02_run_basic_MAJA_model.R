## running MAJA basic models 

eval "$(micromamba shell hook --shell bash )"
micromamba activate BayesR



mpiexec -n 4 python -m mpi4py MAJA/maja.py --n 14537 --p 439 --q 6 --iters 5000 --burnin 1000 --x /GS_MS_prots_MAJA_cog_adPRS_cogPRS_age_sex_kinship_corrected.zarr --y /basic_model/ --diagnostics True --g 439


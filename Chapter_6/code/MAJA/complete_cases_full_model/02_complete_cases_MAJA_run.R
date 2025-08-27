
## complete data  model 

eval "$(micromamba shell hook --shell bash )"
micromamba activate BayesR



mpiexec -n 4 python -m mpi4py MAJA/maja.py --n 11448 --p 439 --q 6 --iters 5000 --burnin 1000 --x /GS_MS_prots_MAJA_cog_adPRS_cogPRS_fully_corrected.zarr --y /GS_cognitive_adPRS_cogPRS_MAJA_full_corrected.txt --dir /complete_data/ --diagnostics True --g 439


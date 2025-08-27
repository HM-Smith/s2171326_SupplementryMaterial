#run imputed covar model in MAJA
eval "$(micromamba shell hook --shell bash )"
micromamba activate MAJA


mpiexec -n 4 python -m mpi4py MAJA/maja.py --n 14258 --p 439 --q 6 --iters 5000 --burnin 1000 --x GS_MS_prots_MAJA_cog_adPRS_cogPRS_KNNIMPUTE_corrected.zarr --y GS_cognitive_adPRS_cogPRS_MAJA_KNNimpute_corrected.txt --dir /KNN_imputed/  --diagnostics True --g 439


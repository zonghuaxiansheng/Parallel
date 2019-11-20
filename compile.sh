mpic++ -g -Wall -std=c++11 src/*.cc -o sim_pipeline -fopenmp -DSIM_PIPELINE
mpic++ -g -Wall -std=c++11 src/*.cc -o sim_alltoall -fopenmp -DSIM_ALLTOALL
mpic++ -g -Wall -std=c++11 src/*.cc -o sim_bcast -fopenmp -DSIM_BCAST

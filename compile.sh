# mpic++ -g -Wall -std=c++11 src/*.cc -o sim_pipeline -fopenmp -DSIM_PIPELINE
# mpic++ -g -Wall -std=c++11 src/*.cc -o sim_alltoall -fopenmp -DSIM_ALLTOALL
# mpic++ -g -Wall -std=c++11 src/*.cc -o sim_bcast -fopenmp -DSIM_BCAST
# mpic++ -g -Wall -std=c++11 src/*.cc -o sim_lusplit -fopenmp -DSIM_LUSPLIT
# mpic++ -g -Wall -std=c++11 src/*.cc -o sim_dishsum -fopenmp -DSIM_DISHSUM
# mpic++ -g -Wall -std=c++11 src/*.cc -o sim_btreesum -fopenmp -DSIM_BTREESUM
# mpic++ -g -Wall -std=c++11 src/*.cc -o sim_pserver -fopenmp -DSIM_PSERVER

# mpic++ -g -Wall -std=c++11 src/*.cc -o sim_mc_single -fopenmp -DSIM_MC_SINGLE
# mpic++ -g -Wall -std=c++11 src/*.cc -o sim_mc_parallel -fopenmp -DSIM_MC_PARALLEL

# mpic++ -g -Wall -std=c++11 src/*.cc -o sim_foxmult -fopenmp -DSIM_FOXMULT

mpic++ -g -Wall -std=c++11 src/*.cc -o sim_tc_parallel -fopenmp -DSIM_TC_PARALLEL

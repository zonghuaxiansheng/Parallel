#include <iostream>
#include "mpi.h"
#include "experiment.h"

int main(int argc, char* argv[]) {

	int my_rank, psize;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &psize);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if (my_rank == 0) {
		std::cout << "Parallel Project by USTC.zonghua !" << std::endl;
	}
#ifdef SIM_PIPELINE
	ustc_parallel::CreatePipeLine(my_rank, psize, MPI_COMM_WORLD);
#endif
#ifdef SIM_ALLTOALL
	ustc_parallel::CreateSimAlltoAll(my_rank, psize, MPI_COMM_WORLD);
#endif
#ifdef SIM_BCAST
	ustc_parallel::CreateSimBcast(my_rank, psize, MPI_COMM_WORLD);
#endif

	MPI_Finalize();

	return 0;
}

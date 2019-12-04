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
	
	// ustc_parallel::CreatePipeLine(my_rank, psize, MPI_COMM_WORLD);
	// ustc_parallel::CreateSimAlltoAll(my_rank, psize, MPI_COMM_WORLD);
	// ustc_parallel::CreateSimBcast(my_rank, psize, MPI_COMM_WORLD);
	// ustc_parallel::CreateLUSplit(my_rank, psize, MPI_COMM_WORLD);

	// ustc_parallel::CreateDishSum(my_rank, psize, MPI_COMM_WORLD);
	// ustc_parallel::CreateBinaryTreeSum(my_rank, psize, MPI_COMM_WORLD);
	// ustc_parallel::CreateParameterServer(my_rank, psize, MPI_COMM_WORLD);
	// ustc_parallel::CreateMonteCarloSingle(my_rank, psize, MPI_COMM_WORLD);
	ustc_parallel::CreateMonteCarloParallel(my_rank, psize, MPI_COMM_WORLD);

	MPI_Finalize();

	return 0;
}
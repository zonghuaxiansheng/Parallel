#ifndef _USTC_EXPERIMENT_H_
#define _USTC_EXPERIMENT_H_

namespace ustc_parallel {

	void CreatePipeLine(int& my_rank, int& psize, MPI_Comm my_comm);
	void CreateSimAlltoAll(int& my_rank, int& psize, MPI_Comm my_comm);
	void CreateSimBcast(int& my_rank, int& psize, MPI_Comm my_comm);
	void CreateLUSplit(int& my_rank, int& psize, MPI_Comm my_comm);

	void CreateDishSum(int& my_rank, int& psize, MPI_Comm my_comm);
	void CreateBinaryTreeSum(int& my_rank, int& psize, MPI_Comm my_comm);
	void CreateFoxMatrixMult(int& my_rank, int& psize, MPI_Comm my_comm);
	void CreateParameterServer(int& my_rank, int& psize, MPI_Comm my_comm);
	void CreateMonteCarloSingle(int& my_rank, int& psize, MPI_Comm my_comm);
	void CreateMonteCarloParallel(int& my_rank, int& psize, MPI_Comm my_comm);
	void CreateTransitiveClosureParallel(int& my_rank, int& psize, MPI_Comm my_comm);
}

#endif // !_USTC_EXPERIMENT_H_

#ifndef _USTC_EXPERIMENT_H_
#define _USTC_EXPERIMENT_H_

namespace ustc_parallel {

	void CreatePipeLine(int& my_rank, int& psize, MPI_Comm my_comm);
	void CreateSimAlltoAll(int& my_rank, int& psize, MPI_Comm my_comm);
	void CreateSimBcast(int& my_rank, int& psize, MPI_Comm my_comm);
}

#endif // !_USTC_EXPERIMENT_H_

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "mpi.h"
#include "utils.h"

namespace ustc_parallel {

double sstart_t = 0.0, send_t = 0.0, pstart_t = 0.0, pend_t = 0.0, trans_t = 0.0;

template<typename T>
void RunTransitiveClosureSerial(T*** A, T*** M, int matrix_size) {
    auto& aptr = *(A);
    auto& mptr = *(M);
    for (int k = 1; k <= log(matrix_size); k ++) {
        for (int i = 0; i < matrix_size; i ++) {
            for (int j = 0; j < matrix_size; j ++) {
                int s = 0;
                while (s < matrix_size && (aptr[i][s] == 0 || aptr[s][j] == 0)) {
                    s ++;
                }
                if (s < matrix_size) {
                    mptr[i][j] = 1;
                } else {
                    mptr[i][j] = 0;
                }
            }
        }
        for (int i = 0; i < matrix_size; i ++) {
            for (int j = 0; j < matrix_size; j ++) {
                aptr[i][j] = mptr[i][j];
            }
        }
        // MatrixPrint<int>(&aptr, "Loop-A", matrix_size);
    }
}

void CreateTransitiveClosureParallel(int& n, int& my_rank, int& psize, MPI_Comm my_comm) {

    int** SA = nullptr;
    int** PA = nullptr;
    int** M = nullptr;
    int sub_matrix_size = n;
    int matrix_size = sub_matrix_size * psize;
    srand(time(NULL));
    CreateMatrix<int>(&SA, matrix_size, matrix_size);
    CreateMatrix<int>(&PA, matrix_size, matrix_size);
    MatrixCopy<int>(&SA, &PA, matrix_size, matrix_size);
    CreateMatrix<int>(&M, matrix_size, matrix_size);

    if (my_rank == 0) {
        sstart_t = MPI_Wtime();
        RunTransitiveClosureSerial<int>(&SA, &M, matrix_size);
        send_t = MPI_Wtime();
    }
    int** LocalA = nullptr;
    int** LocalB = nullptr;
    int** TmpB = nullptr;
    int row_size = matrix_size / psize;
    int col_size = row_size;
    if (matrix_size % psize != 0) {
        std::cout << "* Unsupport matrix_size % psize != 0 now !" << std::endl;
    }

    CreateMatrix<int>(&LocalA, row_size, matrix_size);
    CreateMatrix<int>(&LocalB, matrix_size, col_size);
    CreateMatrix<int>(&TmpB, matrix_size, col_size);

    pstart_t = MPI_Wtime();
    for (int i = 1; i <= log(matrix_size); i ++) {
        // 1. Root send sub_rows and sub_cols to each processor
        if (my_rank == 0) {
            std::vector<MPI_Request> req_vec;
            double tmp_t = MPI_Wtime();
            for (int p = 1; p < psize; p ++) {
                for (int r = 0; r < row_size; r ++) {
                    MPI_Request request;
                    MPI_Isend(PA[p * row_size + r], matrix_size, MPI_INT, p, p, my_comm, &request);
                    req_vec.push_back(request);
                }
                for (int r = 0; r < matrix_size; r ++) {
                    MPI_Request request;
                    MPI_Isend(&PA[r][p * col_size], col_size, MPI_INT, p, p + psize, my_comm, &request);
                    req_vec.push_back(request);
                }
            }
            trans_t += (MPI_Wtime() - tmp_t);
            for (int r = 0; r < row_size; r ++) {
                for (int c = 0; c < matrix_size; c ++) {
                    LocalA[r][c] = PA[r][c];
                }
            }
            for (int r = 0; r < matrix_size; r ++) {
                for (int c = 0; c < col_size; c ++) {
                    LocalB[r][c] = PA[r][c];
                }
            }
            MPI_Status status;
            tmp_t = MPI_Wtime();
            for (auto iter = req_vec.begin(); iter != req_vec.end(); iter ++) {
                MPI_Wait(&(*iter), &status);
            }
            trans_t += (MPI_Wtime() - tmp_t);
        } else {
            MPI_Request request;
            MPI_Status status;
            double tmp_t = MPI_Wtime();
            for (int r = 0; r < row_size; r ++) {
                MPI_Irecv(LocalA[r], matrix_size, MPI_INT, 0, my_rank, my_comm, &request);
                MPI_Wait(&request, &status);
            }
            for (int r = 0; r < matrix_size; r ++) {
                MPI_Irecv(LocalB[r], col_size, MPI_INT, 0, my_rank + psize, my_comm, &request);
                MPI_Wait(&request, &status);
            }
            trans_t += (MPI_Wtime() - tmp_t);
        }
        // 2. Start compute for each processor
        for (int j = 0; j < psize; j ++) {
            for (int r = 0; r < row_size; r ++) {
                for (int c = 0; c < col_size; c ++) {
                    int s = 0;
                    while (s < matrix_size && (LocalA[r][s] == 0 || LocalB[s][c] == 0)) {
                        s ++;
                    }
                    int start_col = (my_rank + j) % psize;
                    M[my_rank * row_size + r][start_col * col_size + c] = (s < matrix_size) ? 1 : 0;
                }
            }
            // Shift LocalB
            std::vector<MPI_Request> req_vec_r;
            std::vector<MPI_Request> req_vec_s;
            MatrixCopy<int>(&LocalB, &TmpB, matrix_size, col_size);
            double tmp_t = MPI_Wtime();
            for (int r = 0; r < matrix_size; r ++) {
                int src = (my_rank + 1) % psize;
                int dest = (my_rank == 0) ? psize - 1 : (my_rank - 1);
                MPI_Request request_s, request_r;
                MPI_Isend(TmpB[r], col_size, MPI_INT, dest, my_rank, my_comm, &request_s);
                MPI_Irecv(LocalB[r], col_size, MPI_INT, src, src, my_comm, &request_r);
                req_vec_s.push_back(request_s);
                req_vec_r.push_back(request_r);	
            }
            MPI_Status status;
            for (auto iter = req_vec_s.begin(); iter != req_vec_s.end(); iter ++) {
                MPI_Wait(&(*iter), &status);
            }
            for (auto iter = req_vec_r.begin(); iter != req_vec_r.end(); iter ++) {
                MPI_Wait(&(*iter), &status);
            }
            trans_t += (MPI_Wtime() - tmp_t);
        }
        // 3. Collect compute output
        if (my_rank != 0) {
            // Send sub_rows to root
            double tmp_t = MPI_Wtime();
            for (int r = 0; r < row_size; r ++) {
                MPI_Send(M[my_rank * row_size + r], matrix_size, MPI_INT, 0, my_rank, my_comm);
            }
            trans_t += (MPI_Wtime() - tmp_t);
        } else {
            // Recv sub_rows from other processors
            double tmp_t = MPI_Wtime();
            for (int p = 1; p < psize; p ++) {
                for (int r = 0; r < row_size; r ++) {
                    MPI_Status status;
                    MPI_Recv(M[p * row_size + r], matrix_size, MPI_INT, p, p, my_comm, &status);
                }	
            }
            trans_t += (MPI_Wtime() - tmp_t);
        }
        // 4. Copy M to A
        if (my_rank == 0) {
            for (int i = 0; i < matrix_size; i ++) {
                for (int j = 0; j < matrix_size; j ++) {
                    PA[i][j] = M[i][j];
                }
            }
        }
    }
    pend_t = MPI_Wtime();
    // Root show matrix A
    if (my_rank == 0) {
        // MatrixPrint<int>(&PA, "Parallel-A", matrix_size);
        MatrixCompare<int>(&SA, &PA, matrix_size);
        std::cerr << n << "\t" << psize << "\t" << (send_t - sstart_t) << "\t" 
            << (pend_t - pstart_t) << "\t" << trans_t << "\t" << (pend_t - pstart_t - trans_t) << std::endl;
    }

    ReleaseMatrix<int>(&SA, matrix_size);
    ReleaseMatrix<int>(&PA, matrix_size);
    ReleaseMatrix<int>(&M, matrix_size);
    ReleaseMatrix<int>(&LocalA, row_size);
    ReleaseMatrix<int>(&LocalB, matrix_size);
    ReleaseMatrix<int>(&TmpB, matrix_size);
}

}

int main(int argc, char* argv[]) {

    int my_rank, psize;
    int n = 4;
    
    if (argc == 2) {
        n = atoi(argv[1]);
    }

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &psize);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);   

    ustc_parallel::CreateTransitiveClosureParallel(n, my_rank, psize, MPI_COMM_WORLD);

    MPI_Finalize();

	return 0;
}
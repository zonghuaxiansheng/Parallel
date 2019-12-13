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
void RunGaussEliminSerial(T*** A, T** B, T** X, int row_size, int col_size) {
    auto& aptr = *(A);
    auto& bptr = *(B);
    auto& xptr = *(X);
    int* shift = new int[col_size];
    for (int i = 0; i < col_size; i ++) {
        shift[i] = i;
    }
    for (int k = 0; k < row_size; k ++) {
        T max_elem = 0;
        int elem_row = 0;
        int elem_col = 0;
        // Select max item
        for (int r = k; r < row_size; r ++) {
            for (int c = k; c < col_size; c ++) {
                if (fabs(aptr[r][c]) > max_elem) {
                    max_elem = fabs(aptr[r][c]);
                    elem_row = r;
                    elem_col = c;
                }
            }
        }
        // Switch max col
        if (elem_col != k) {
            for (int r = 0; r < row_size; r ++) {
                auto buf = aptr[r][k];
                aptr[r][k] = aptr[r][elem_col];
                aptr[r][elem_col] = buf;
            }
            auto buf = shift[k];
            shift[k] = shift[elem_col];
            shift[elem_col] = buf;
        }
        // Switch max row
        if (elem_row != k) {
            for (int c = 0; c < col_size; c ++) {
                auto buf = aptr[k][c];
                aptr[k][c] = aptr[elem_row][c];
                aptr[elem_row][c] = buf;
            }
            auto buf = bptr[k];
            bptr[k] = bptr[elem_row];
            bptr[elem_row] = buf;
        }
        // Div current line
        for (int c = k + 1; c < col_size; c ++) {
            aptr[k][c] = aptr[k][c] / aptr[k][k];
        }
        bptr[k] = bptr[k] / aptr[k][k];
        aptr[k][k] = 1;
        // Start eliminate
        for (int r = k + 1; r < row_size; r ++) {
            for (int c = k + 1; c < col_size; c ++) {
                aptr[r][c] -= aptr[r][k] * aptr[k][c];
            }
            bptr[r] -= aptr[r][k] * bptr[k];
            aptr[r][k] = 0;
        }
    }
    for (int c = col_size - 1; c >= 0; c --) {
        if (c == col_size - 1) {
            xptr[c] = bptr[c] / aptr[c][c];
        } else {
            float sub_sum = 0;
            for (int c1 = c + 1; c1 < col_size; c1 ++) {
                sub_sum += aptr[c][c1] * xptr[c1];
            }
            xptr[c] = (bptr[c] - sub_sum) / aptr[c][c];
        }
    }
    T* SortX = nullptr;
    CreateVector<T>(&SortX, col_size);
    for (int r = 0; r < col_size; r ++) {
        for (int i = 0; i < col_size; i ++) {
            if (shift[i] == r) {
                SortX[r] = xptr[i];
            }
        }
    }
    VectorCopy<T>(&SortX, &xptr, col_size);
    ReleaseVector<T>(&SortX);
}

void CreateGaussEliminParallel(int& n, int& my_rank, int& psize, MPI_Comm my_comm) {
    float** SA = nullptr;
    float** PA = nullptr;
    float* SB = nullptr;
    float* PB = nullptr;
    float* SX = nullptr;
    float* PX = nullptr;
    int sub_matrix_size = n;
    int row_size = sub_matrix_size * psize;
    int col_size = sub_matrix_size * psize;
    CreateMatrix<float>(&PA, row_size, col_size, init_t::GS);
    CreateVector<float>(&PB, row_size, init_t::GS);
    CreateVector<float>(&PX, col_size);
    if (my_rank == 0) {
        CreateMatrix<float>(&SA, row_size, col_size, init_t::GS);
        CreateVector<float>(&SB, row_size, init_t::GS);
        CreateVector<float>(&SX, col_size);
        MatrixCopy<float>(&PA, &SA, row_size, col_size);
        VectorCopy<float>(&PB, &SB, row_size);
        sstart_t = MPI_Wtime();
        RunGaussEliminSerial<float>(&SA, &SB, &SX, row_size, col_size);
        send_t = MPI_Wtime();
        ReleaseMatrix<float>(&SA, row_size);
        ReleaseVector<float>(&SB);
    }
    MPI_Barrier(my_comm);

    float** LocalA = nullptr;
    float* LocalB = nullptr;
    float* GlobalF = nullptr;
    CreateMatrix<float>(&LocalA, sub_matrix_size, col_size);
    CreateVector<float>(&LocalB, sub_matrix_size);
    CreateVector<float>(&GlobalF, col_size + 1);

    pstart_t = MPI_Wtime();
    for (int r = 0; r < sub_matrix_size; r ++) {
        VectorCopy<float>(&PA[r * psize + my_rank], &LocalA[r], col_size);
        LocalB[r] = PB[r * psize + my_rank];
    }
    int* shift = new int[col_size];
    for (int i = 0; i < col_size; i ++) {
        shift[i] = i;
    }

    float lmax[4] = {0};
    float* lmax_all_put = new float[psize * 4];
    float* lmax_all_get = new float[psize * 4];
    // 1. forward
    for (int s = 0; s < sub_matrix_size; s ++) {
        for (int p = 0; p < psize; p ++) {
            // 1.1 find local max item
            int mainIdx = s * psize + p;
            if (my_rank < p) {
                lmax[0] = 0;
                for (int ss = s + 1; ss < sub_matrix_size; ss ++) {
                    for (int sc = mainIdx; sc < col_size; sc ++) {
                        if (fabs(LocalA[ss][sc]) > lmax[0]) {
                            lmax[0] = fabs(LocalA[ss][sc]);
                            lmax[1] = ss;
                            lmax[2] = sc;
                            lmax[3] = my_rank;
                        }
                    }
                }
            } else {
                lmax[0] = 0;
                for (int ss = s; ss < sub_matrix_size; ss ++) {
                    for (int sc = mainIdx; sc < col_size; sc ++) {
                        if (fabs(LocalA[ss][sc]) > lmax[0]) {
                            lmax[0] = fabs(LocalA[ss][sc]);
                            lmax[1] = ss;
                            lmax[2] = sc;
                            lmax[3] = my_rank;
                        }
                    }
                }
            }
            // 1.3 copy local max item
            for (int i = 0; i < psize; i ++) {
                for (int j = 0; j < 4; j ++) {
                    lmax_all_put[i * 4 + j] = lmax[j];
                }
            }
            double tmp_t = MPI_Wtime();
            MPI_Alltoall(lmax_all_put, 4, MPI_FLOAT, lmax_all_get, 4, MPI_FLOAT, my_comm);
            trans_t += (MPI_Wtime() - tmp_t);
            // 1.4 compute global max item
            lmax[0] = 0;
            for (int i = 0; i < psize; i ++) {
                if (fabs(lmax_all_get[i * 4]) > fabs(lmax[0])) {
                    for (int j = 0; j < 4; j ++) {
                        lmax[j] = lmax_all_get[i * 4 + j];
                    }
                }
            }
            // 1.5 switch col
            if ((int)lmax[2] != mainIdx) {
                for (int r = 0; r < sub_matrix_size; r ++) {
                    auto buf = LocalA[r][mainIdx];
                    LocalA[r][mainIdx] = LocalA[r][(int)lmax[2]];
                    LocalA[r][(int)lmax[2]] = buf;
                }
                auto buf = shift[mainIdx];
                shift[mainIdx] = shift[(int)lmax[2]];
                shift[(int)lmax[2]] = buf;
            }
            // 1.6 switch row
            if ((my_rank == p) && (lmax[3] == p)) {
                if ((int)lmax[1] != s) {
                    for (int c = 0; c < col_size; c ++) {
                        auto buf = LocalA[s][c];
                        LocalA[s][c] = LocalA[(int)lmax[1]][c];
                        LocalA[(int)lmax[1]][c] = buf;
                    }
                    auto buf = LocalB[s];
                    LocalB[s] = LocalB[(int)lmax[1]];
                    LocalB[(int)lmax[1]] = buf;
                }
            } else if (my_rank == p || my_rank == (int)lmax[3]) {
                float* putTmpA = nullptr;
                float* getTmpA = nullptr;
                CreateVector<float>(&putTmpA, col_size + 1);
                CreateVector<float>(&getTmpA, col_size + 1);
                MPI_Status status;
                if (my_rank == p) {
                    // recv
                    VectorCopy<float>(&LocalA[s], &putTmpA, col_size);
                    putTmpA[col_size] = LocalB[s];
                    tmp_t = MPI_Wtime();
                    MPI_Recv(getTmpA, col_size + 1, MPI_FLOAT, (int)lmax[3], (int)lmax[3], my_comm, &status);
                    MPI_Send(putTmpA, col_size + 1, MPI_FLOAT, (int)lmax[3], my_rank, my_comm);
                    trans_t += (MPI_Wtime() - tmp_t);
                    VectorCopy<float>(&getTmpA, &LocalA[s], col_size);
                    LocalB[s] = getTmpA[col_size];
                }
                if (my_rank == (int)lmax[3]) {
                    // send
                    VectorCopy<float>(&LocalA[(int)lmax[1]], &putTmpA, col_size);
                    putTmpA[col_size] = LocalB[(int)lmax[1]];
                    tmp_t = MPI_Wtime();
                    MPI_Send(putTmpA, col_size + 1, MPI_FLOAT, p, my_rank, my_comm);
                    MPI_Recv(getTmpA, col_size + 1, MPI_FLOAT, p, p, my_comm, &status);
                    trans_t += (MPI_Wtime() - tmp_t);
                    VectorCopy<float>(&getTmpA, &LocalA[(int)lmax[1]], col_size);
                    LocalB[(int)lmax[1]] = getTmpA[col_size];
                }
                ReleaseVector<float>(&putTmpA);
                ReleaseVector<float>(&getTmpA);						
            }
            // 1.7 compute main line
            if (my_rank == p) {
                for (int c = mainIdx + 1; c < col_size; c ++) {
                    LocalA[s][c] = LocalA[s][c] / LocalA[s][mainIdx];
                }
                LocalB[s] = LocalB[s] / LocalA[s][mainIdx];
                LocalA[s][mainIdx] = 1;
                VectorCopy<float>(&LocalA[s], &GlobalF, col_size);
                GlobalF[col_size] = LocalB[s];
            }
            tmp_t = MPI_Wtime();
            MPI_Bcast(GlobalF, col_size + 1, MPI_FLOAT, p, my_comm);
            trans_t += (MPI_Wtime() - tmp_t);
            // 1.8 compute sub matrix
            int start_row = (my_rank <= p) ? s + 1 : s;
            for (int r = start_row; r < sub_matrix_size; r ++) {
                for (int c = mainIdx + 1; c < col_size; c ++) {
                    LocalA[r][c] -= GlobalF[c] * LocalA[r][mainIdx];
                }
                LocalB[r] -= GlobalF[col_size] * LocalA[r][mainIdx];
                LocalA[r][mainIdx] = 0;
            }
        }
    }
    MPI_Barrier(my_comm);
    // 2. backward
    float* sum = new float[sub_matrix_size];
    for (int s = 0; s < sub_matrix_size; s ++) {
        sum[s] = 0.0;
    }
    for (int s = sub_matrix_size - 1; s >= 0; s --) {
        for (int p = psize - 1; p >= 0; p --) {
            if (my_rank == p) {
                PX[s * psize + p] = (LocalB[s] - sum[s]) / LocalA[s][s * psize + p];
                double tmp_t = MPI_Wtime();
                MPI_Bcast(&PX[s * psize + p], 1, MPI_FLOAT, p, my_comm);
                trans_t += (MPI_Wtime() - tmp_t);
                for (int ss = 0; ss < sub_matrix_size; ss ++) {
                    sum[ss] += LocalA[ss][s * psize + p] * PX[s * psize + p];
                }
            } else {
                double tmp_t = MPI_Wtime();
                MPI_Bcast(&PX[s * psize + p], 1, MPI_FLOAT, p, my_comm);
                trans_t += (MPI_Wtime() - tmp_t);
                int end_row = (my_rank > p) ? s - 1 : s;
                for (int ss = 0; ss <= end_row; ss ++) {
                    sum[ss] += LocalA[ss][s * psize + p] * PX[s * psize + p];
                }
            }
            MPI_Barrier(my_comm);
        }
    }
    pend_t = MPI_Wtime();
    if (my_rank == 0) {
        float* SortX = nullptr;
        CreateVector<float>(&SortX, col_size);
        for (int r = 0; r < col_size; r ++) {
            for (int i = 0; i < col_size; i ++) {
                if (shift[i] == r) {
                    SortX[r] = PX[i];
                }
            }
        }
        VectorCopy<float>(&SortX, &PX, col_size);
        if (VectorCompare<float>(&SX, &PX, col_size)) {
            std::cout << "* Serial vs Parallel pass !" << std::endl;
        }
        std::cerr << n << "\t" << psize << "\t" << (send_t - sstart_t) << "\t" 
            << (pend_t - pstart_t) << "\t" << trans_t << "\t" << (pend_t - pstart_t - trans_t) << std::endl;
        ReleaseVector<float>(&SX);
        ReleaseVector<float>(&SortX);
    }

    delete [] shift;
    delete [] lmax_all_get;
    delete [] lmax_all_put;
    delete [] sum;
    ReleaseMatrix<float>(&PA, row_size);
    ReleaseMatrix<float>(&LocalA, sub_matrix_size);
    ReleaseVector<float>(&PB);
    ReleaseVector<float>(&LocalB);
    ReleaseVector<float>(&PX);
    ReleaseVector<float>(&GlobalF);
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

    ustc_parallel::CreateGaussEliminParallel(n, my_rank, psize, MPI_COMM_WORLD);

    MPI_Finalize();

	return 0;
}

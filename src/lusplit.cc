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
#include "omp.h"
#include "utils.h"

namespace ustc_parallel {

double sstart_t = 0.0, send_t = 0.0, pstart_t = 0.0, pend_t = 0.0, trans_t = 0.0;

template<typename T>
void DumpMatrix(T** *mptr_, std::string mname, int matrix_size) {
    std::ofstream ofile;
    ofile.open(mname, std::ios::out);
    if (!ofile.is_open()) {
        std::cout << "* Open file failed!" << std::endl;
        exit(1);
    }
    auto& mptr = *(mptr_);
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            ofile << mptr[i][j] << " ";
        }
    }
    ofile.close();
}

template<typename T>
void LoadMatrix(T** *mptr_, std::string mname, int matrix_size) {
    std::ifstream ifile;
    ifile.open(mname, std::ios::in);
    if (!ifile.is_open()) {
        std::cerr << "* Open file failed!" << std::endl;
        exit(1);
    } else {
        std::cout << "* Open file " << mname << std::endl;
    }
    auto& mptr = *(mptr_);
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            ifile >> mptr[i][j];
        }
    }
    ifile.close();
}

template<typename T>
void PrintMatrix(T** mptr, std::string name, int matrix_size) {
    std::cout << std::endl << "* Matrix " << name << " ..." << std::endl;
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            std::cout << std::setw(5) << mptr[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

template<typename T>
void CreateLuMatrix(T** *L, T** *U, T** *A, int matrix_size) {
    auto& lptr = *(L);
    auto& uptr = *(U);
    auto& luptr = *(A);
    lptr = new T * [matrix_size];
    uptr = new T * [matrix_size];
    luptr = new T * [matrix_size];

    for (int i = 0; i < matrix_size; i++) {
        lptr[i] = new T[matrix_size];
        uptr[i] = new T[matrix_size];
        luptr[i] = new T[matrix_size];
        memset(luptr[i], 0, matrix_size * sizeof(T));
    }

    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            if (i < j) {
                uptr[i][j] = rand() % 50 + 1;
                lptr[i][j] = 0;
            } else if (i == j) {
                uptr[i][j] = rand() % 50 + 1;
                lptr[i][j] = 1;
            } else {
                uptr[i][j] = 0;
                lptr[i][j] = rand() % 50 + 1;
            }
        }			
    }

    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            for (int k = 0; k < matrix_size; k++) {
                luptr[i][j] += lptr[i][k] * uptr[k][j];
            }
        }

    }
}

template<typename T>
void SlpitA2LU(T** *A, T** *L, T** *U, int matrix_size) {
    auto& lptr = *(L);
    auto& uptr = *(U);
    auto& luptr = *(A);

    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            if (i < j) {
                uptr[i][j] = luptr[i][j];
                lptr[i][j] = 0;
            }
            else if (i == j) {
                uptr[i][j] = luptr[i][j];
                lptr[i][j] = 1;
            }
            else {
                uptr[i][j] = 0;
                lptr[i][j] = luptr[i][j];
            }
        }
    }
}

template<typename T>
void ReleaseLuMatrix(T** lptr, T** uptr, T** luptr, int matrix_size) {
    for (int i = 0; i < matrix_size; i++) {
        delete[] lptr[i];
        delete[] uptr[i];
        delete[] luptr[i];
    }
    delete[] lptr;
    delete[] uptr;
    delete[] luptr;
}

void CreateLUSplit(int& my_rank, int& psize, MPI_Comm my_comm) {

    int **A = nullptr, **L = nullptr, **U = nullptr;
    int matrix_size = 500;
    CreateLuMatrix<int>(&L, &U, &A, matrix_size);

    // if (my_rank == 0) {
    // 	DumpMatrix<int>(&L, "Matrix_L", matrix_size);
    // 	DumpMatrix<int>(&U, "Matrix_U", matrix_size);
    // 	DumpMatrix<int>(&A, "Matrix_A", matrix_size);
    // 	PrintMatrix<int>(L, "L", matrix_size);
    // 	PrintMatrix<int>(U, "U", matrix_size);
    // 	PrintMatrix<int>(A, "A", matrix_size);
    // }
    // MPI_Barrier(my_comm);
    for (int i = 0; i < psize; i++) {
        if (my_rank == i) {
            LoadMatrix<int>(&L, "Matrix_L", matrix_size);
            LoadMatrix<int>(&U, "Matrix_U", matrix_size);
            LoadMatrix<int>(&A, "Matrix_A", matrix_size);
        }
        MPI_Barrier(my_comm);
    }

    int turn = ceil((float)(matrix_size) / (float)(psize));
    int* mainPtr = new int[matrix_size];
    memset(mainPtr, 0, matrix_size * sizeof(int));

    pstart_t = MPI_Wtime();
    // For each trun
    for (int t = 0; t < turn; t++) {
        // For each processor
        for (int p = 0; p < psize; p++) {
            // My trun
            auto mainLine = t * psize + p;
            if (mainLine < matrix_size) {
                if (my_rank == p) {
                    // std::cout << std::endl << "* Main Line " << mainLine << ": ";
                    for (int k = 0; k < matrix_size; k++) {
                        mainPtr[k] = A[mainLine][k];
                        // std::cout << mainPtr[k] << " ";
                    }
                    // std::cout << std::endl;
                }
                double tmp_t = MPI_Wtime();
                MPI_Bcast(mainPtr, matrix_size, MPI_INT, p, my_comm);
                MPI_Barrier(my_comm);
                trans_t += (MPI_Wtime() - tmp_t);
                for (int j = t; j < turn; j++) {
                    int curLine = j * psize + my_rank;
                    if (curLine > mainLine && curLine < matrix_size) {
                        A[curLine][mainLine] /= mainPtr[mainLine];
                        for (int i = mainLine + 1; i < matrix_size; i++) {
                            A[curLine][i] -= A[curLine][mainLine] * mainPtr[i];
                        }
                    }
                }
                if (my_rank != p) {
                    for (int i = 0; i < matrix_size; i++) {
                        A[mainLine][i] = mainPtr[i];
                    }
                }
            }
        }
    }
    // MPI_Barrier(my_comm);
    pend_t = MPI_Wtime();
    if (my_rank == 0) {
        // PrintMatrix<int>(A, "Final A", matrix_size);
        int** AS = nullptr, ** LS = nullptr, ** US = nullptr;
        CreateLuMatrix<int>(&LS, &US, &AS, matrix_size);
        SlpitA2LU<int>(&A, &LS, &US, matrix_size);

        for (int i = 0; i < matrix_size; i++) {
            for (int j = 0; j < matrix_size; j++) {
                if (L[i][j] != LS[i][j]) {
                    std::cerr << "LU Split error with L<" << i << "," << j << ">(" 
                            << L[i][j] << "," << LS[i][j] << ")" << std::endl;
                    exit(1);
                }
                if (U[i][j] != US[i][j]) {
                    std::cerr << "LU Split error with U<" << i << "," << j << ">("
                        << U[i][j] << "," << US[i][j] << ")" << std::endl;
                    exit(1);
                }
            }
        }
        
        std::cerr << psize << "\t" << (pend_t - pstart_t) << "\t" 
            << trans_t << "\t" << (pend_t - pstart_t - trans_t) << std::endl;

        // PrintMatrix<int>(LS, "Final L", matrix_size);
        // PrintMatrix<int>(US, "Final U", matrix_size);
    }

    // ReleaseLuMatrix<int>(L, U, A, matrix_size);
}

}

int main(int argc, char* argv[]) {

    int my_rank, psize;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &psize);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);   

    ustc_parallel::CreateLUSplit(my_rank, psize, MPI_COMM_WORLD);

    MPI_Finalize();

	return 0;
}

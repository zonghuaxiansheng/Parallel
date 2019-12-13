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
#include "experiment.h"
#ifdef WINDOWS_PLATFORM
#include <Windows.h>
#else
#include <unistd.h>
#endif	// WINDOWS_PLATFORM

namespace ustc_parallel {

	struct Node
	{
		int key_;
		int color_;
	};

	void inline gethostname(int my_rank, char* hostname, int len) {
		char tname[16];
		int group = my_rank % 3;
		if (group == 0) {
			strcpy(tname, "Node0");
		} else if (group == 1) {
			strcpy(tname, "Node1");
		} else {
			strcpy(tname, "Node2");
		}
		strcpy(hostname, tname);
	}

	void CreatePipeLine(int& my_rank, int& psize, MPI_Comm my_comm) {
		MPI_Request request_r, request_s;
		MPI_Status status_r;
		int pre_rank = my_rank - 1;
		int next_rank = my_rank + 1;
		int tag = 10;
		char next_rbuf[128];
		char cur_rbuf[128];
		char pre_sbuf[128] = "A";
		char cur_sbuf[128] = "A";
		int buf_len = 1;
		int rcnt = 0;

		while (rcnt < 10) {
			if (my_rank == 0) {
				// MPI_Send(sbuf, 128, MPI_CHAR, next_rank, tag, my_comm);
				MPI_Isend(pre_sbuf, 128, MPI_CHAR, next_rank, tag, my_comm, &request_s);
				// MPI_Wait(&request_s, &status_s);
				std::cout << "Rank: " << my_rank
					<< " Sbuf: " << pre_sbuf << " Dst Rank: " << next_rank << std::endl;
			} else if (my_rank == psize - 1) {
				MPI_Irecv(next_rbuf, 128, MPI_CHAR, pre_rank, tag, my_comm, &request_r);
				MPI_Wait(&request_r, &status_r);
				std::cout << "Rank: " << my_rank
					<< " Rbuf: " << next_rbuf << " Src Rank: " << status_r.MPI_SOURCE << std::endl;
			} else {
				MPI_Irecv(next_rbuf, 128, MPI_CHAR, pre_rank, tag, my_comm, &request_r);
				MPI_Isend(pre_sbuf, 128, MPI_CHAR, next_rank, tag, my_comm, &request_s);
				MPI_Wait(&request_r, &status_r);
				// MPI_Wait(&request_s, &status_s);
				std::cout << "Rank: " << my_rank
					<< " Rbuf: " << next_rbuf << " Src Rank: " << status_r.MPI_SOURCE << std::endl;
				strcpy(cur_rbuf, next_rbuf);
				for (int i = 0; i < buf_len; i++) {
					cur_sbuf[i] = cur_rbuf[i] + 1;
				}
				cur_sbuf[buf_len] = '\0';
				strcpy(pre_sbuf, cur_sbuf);
				// MPI_Send(sbuf, 128, MPI_CHAR, next_rank, tag, my_comm);
				std::cout << "Rank: " << my_rank
					<< " Sbuf: " << pre_sbuf << " Dst Rank: " << next_rank << std::endl;
			}
			rcnt ++;
		}
	}

	void CreateSimAlltoAll(int& my_rank, int& psize, MPI_Comm my_comm) {
		MPI_Status status;
		int tag = 10;
		char sbuf[16] = "abc";
		char* rbuf = new char[psize*16];
		memset(rbuf, '\0', psize * 16);

		for (int i = 0; i < psize; i++) {
			if (my_rank != i) {
				MPI_Send(sbuf, 16, MPI_CHAR, i, tag, my_comm);
			}
		}

		for (int i = 0; i < psize; i++) {
			if (my_rank != i) {
				MPI_Recv(rbuf + i * 16, 16, MPI_CHAR, i, tag, my_comm, &status);
			}
		}

		strcpy(rbuf + my_rank * 16, sbuf);

		for (int i = 0; i < psize; i++) {
			char tbuf[16];
			strcpy(tbuf, rbuf + i * 16);
			std::cout << "Rank: " << my_rank
				<< " Buf[" << i << "]: " << tbuf << std::endl;
		}

		delete [] rbuf;
	}

	void CreateSimBcast(int& my_rank, int& psize, MPI_Comm my_comm) {
	
		char hostname[16];
		gethostname(my_rank, hostname, 16);

		char* allname = new char[psize * 16];
		Node* nodes = new Node[psize];
		memset(allname, '\0', psize * 16);

		MPI_Gather(hostname, 16, MPI_CHAR, allname, 16, MPI_CHAR, 0, MPI_COMM_WORLD);

		int sub_color = 0, sub_key = 0;
		int color_tag = 10, key_tag = 20;

		if (my_rank == 0) {
			std::map<std::string, std::pair<int, int>> name_map;
			int sub_color = 0;
			for (int i = 0; i < psize; i++) {
				char tname[16];
				strcpy(tname, allname + i * 16);
				std::string name_str(tname);
				// std::cout << tname << std::endl;
				std::map<std::string, std::pair<int, int>>::iterator iter = name_map.find(name_str);
				if (iter != name_map.end()) {
					nodes[i].color_ = iter->second.first;
					nodes[i].key_ = ++ iter->second.second;
				}
				else {
					nodes[i].color_ = sub_color;
					nodes[i].key_ = 0;
					name_map.insert(std::make_pair(name_str, std::make_pair(sub_color++, 0)));
				}
			}
			for (int i = 0; i < psize; i++) {
				// std::cout << "Node: " << i << " " << nodes[i].color_ << " " << nodes[i].key_ << std::endl;
				if (my_rank != i) {
					MPI_Send(&nodes[i].color_, 1, MPI_INT, i, color_tag, my_comm);
					MPI_Send(&nodes[i].key_, 1, MPI_INT, i, key_tag, my_comm);
				}
			}
			sub_color = nodes[0].color_;
			sub_key = nodes[0].key_;
		} else {
			MPI_Status status;
			MPI_Recv(&sub_color, 1, MPI_INT, 0, color_tag, my_comm, &status);
			MPI_Recv(&sub_key, 1, MPI_INT, 0, key_tag, my_comm, &status);
		}

		std::cout << "Rank: " << my_rank << " Sub Color: " << sub_color << " Sub Key: " << sub_key << std::endl;

		MPI_Barrier(my_comm);

		MPI_Comm leader_comm = MPI_COMM_WORLD, inside_comm = MPI_COMM_WORLD;
		MPI_Comm_split(my_comm, sub_color, sub_key, &inside_comm);
		MPI_Comm_split(my_comm, sub_key, sub_color, &leader_comm);
		
		char sbuf[16];
		if (my_rank == 0) {
			strcpy(sbuf, "sim_bcast");
		} else {
			strcpy(sbuf, "");
		}
		
		std::cout << "Before Split -> Rank: " << my_rank << " Sbuf: " << sbuf << std::endl;
		MPI_Barrier(my_comm);

		MPI_Bcast(sbuf, 16, MPI_CHAR, 0, leader_comm);
		std::cout << "Leader Split -> Rank: " << my_rank << " Sbuf: " << sbuf << std::endl;
		MPI_Barrier(my_comm);

		MPI_Bcast(sbuf, 16, MPI_CHAR, 0, inside_comm);
		std::cout << "Inside Split -> Rank: " << my_rank << " Sbuf: " << sbuf << std::endl;
	}

	// 4. LU Split
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

		// For each trun
		for (int t = 0; t < turn; t++) {
			// For each processor
			for (int p = 0; p < psize; p++) {
				// My trun
				auto mainLine = t * psize + p;
				if (mainLine < matrix_size) {
					if (my_rank == p) {
						std::cout << std::endl << "* Main Line " << mainLine << ": ";
						for (int k = 0; k < matrix_size; k++) {
							mainPtr[k] = A[mainLine][k];
							std::cout << mainPtr[k] << " ";
						}
						std::cout << std::endl;
					}
					MPI_Bcast(mainPtr, matrix_size, MPI_INT, p, my_comm);
					MPI_Barrier(my_comm);
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
		MPI_Barrier(my_comm);

		if (my_rank == 0) {
			PrintMatrix<int>(A, "Final A", matrix_size);
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
			
			PrintMatrix<int>(LS, "Final L", matrix_size);
			PrintMatrix<int>(US, "Final U", matrix_size);
		}

		// ReleaseLuMatrix<int>(L, U, A, matrix_size);
	}

	// 1. Sum
	// 1.1. Dish
	void CreateDishSum(int& my_rank, int& psize, MPI_Comm my_comm) {
		
		auto& n = psize;
		int* sumArray = new int[n];

		int check_sum = 0;
		for (int i = 0; i < n; i++) {
			sumArray[i] = i;
			check_sum += i;
		}

		clock_t start_t = clock();
		// Step
		for (int s = 0; pow(2, s + 1) <= n; s ++) {
			int barrier = pow(2, s);
			int group = pow(2, s + 1);

			auto index = my_rank % group;
			if (index < barrier) {
				int dest = my_rank + barrier;
				MPI_Send(&sumArray[my_rank], 1, MPI_INT, dest, s, my_comm);
			}
			else {
				int dest = my_rank - barrier;
				MPI_Send(&sumArray[my_rank], 1, MPI_INT, dest, s, my_comm);
			}
			// MPI_Barrier(my_comm);
			MPI_Status status;
			if (index < barrier) {
				int dest = my_rank + barrier;
				MPI_Recv(&sumArray[dest], 1, MPI_INT, dest, s, my_comm, &status);
				sumArray[my_rank] += sumArray[dest];
			}
			else {
				int dest = my_rank - barrier;
				MPI_Recv(&sumArray[dest], 1, MPI_INT, dest, s, my_comm, &status);
				sumArray[my_rank] += sumArray[dest];
			}
			// MPI_Barrier(my_comm);
		}
		clock_t end_t = clock();
		int use_t = (end_t - start_t);
		int* use_time_arr = new int[psize];
		MPI_Gather(&use_t, 1, MPI_INT, use_time_arr, 1, MPI_INT, 0, my_comm);

		if (my_rank == 0) {
			float mpi_sum_t = 0.0;
			for (int i = 0; i < psize; i ++) {
				mpi_sum_t += use_time_arr[i];
			}
			mpi_sum_t = mpi_sum_t / (float)psize;
			std::cerr << n << "\t" << mpi_sum_t << std::endl;
		}

		for (int i = 0; i < n; i++) {
			if (my_rank == i ) {
				if (sumArray[i] != check_sum) {
					std::cout << "* Rank: " << my_rank << " Sum Check failed !" << std::endl;
				}
			}
		}
	}
	// 1.2. Binary Tree
	void CreateBinaryTreeSum(int& my_rank, int& psize, MPI_Comm my_comm) {

		auto& n = psize;
		int* sumArray = new int[n];

		int check_sum = 0;
		for (int i = 0; i < n; i++) {
			sumArray[i] = i;
			check_sum += i;
		}
		clock_t start_t = clock();
		// Step
		for (int s = 0; pow(2, s + 1) <= n; s++) {
			int div = pow(2, s + 1);
			int sub = pow(2, s);
			
			if ((my_rank % sub) == 0) {
				if ((my_rank % div) != 0) {
					int dest = my_rank - sub;
					MPI_Send(&sumArray[my_rank], 1, MPI_INT, dest, s, my_comm);
				} else {
					int src = my_rank + sub;
					MPI_Status status;
					MPI_Recv(&sumArray[src], 1, MPI_INT, src, s, my_comm, &status);
					sumArray[my_rank] += sumArray[src];
				}
			}
		}
		MPI_Bcast(&sumArray[0], 1, MPI_INT, 0, my_comm);
		clock_t end_t = clock();
		int use_t = (end_t - start_t);
		int* use_time_arr = new int[psize];
		MPI_Gather(&use_t, 1, MPI_INT, use_time_arr, 1, MPI_INT, 0, my_comm);

		if (my_rank == 0) {
			float mpi_sum_t = 0.0;
			for (int i = 0; i < psize; i ++) {
				mpi_sum_t += use_time_arr[i];
			}
			mpi_sum_t = mpi_sum_t / (float)psize;
			std::cerr << n << "\t" << mpi_sum_t << std::endl;
		}

		if (sumArray[0] != check_sum) {
			std::cout << "* Rank: " << my_rank << " Sum Check failed !" << std::endl;
		}
	}

	template<typename T>
	void MatrixMult(T*** X, T*** Y,T*** OUT, int matrix_size) {
		auto& xptr = *(X);
		auto& yptr = *(Y);
		auto& optr = *(OUT);
		for (int i = 0; i < matrix_size; i ++) {
			for (int j = 0; j < matrix_size; j ++) {
				optr[i][j] = 0;
				for (int k = 0; k < matrix_size; k ++) {
					optr[i][j] += xptr[i][k] * yptr[k][j];
				}
			}
		}
	}

	template<typename T>
	void MatrixAdd(T*** X, T*** Y,T*** OUT, int matrix_size) {
		auto& xptr = *(X);
		auto& yptr = *(Y);
		auto& optr = *(OUT);
		for (int i = 0; i < matrix_size; i ++) {
			for (int j = 0; j < matrix_size; j ++) {
				optr[i][j] = xptr[i][j] + yptr[i][j];
			}
		}
	}

	template<typename T>
	void MatrixPrint(T*** X, std::string name, int row_size, int col_size=0) {
		auto& xptr = *(X);
		std::cout << "* MatrixPrint " << name << std::endl;
		col_size = (col_size == 0) ? row_size : col_size;
		for (int i = 0; i < row_size; i ++) {
			for (int j = 0; j < col_size; j ++) {
				std::cout << xptr[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	template<typename T>
	void MatrixCompare(T*** X, T*** Y, int matrix_size) {
		auto& xptr = *(X);
		auto& yptr = *(Y);
		for (int i = 0; i < matrix_size; i ++) {
			for (int j = 0; j < matrix_size; j ++) {
				if (xptr[i][j] != yptr[i][j]) {
					std::cout << "* Matrix compare failed with " 
							  << xptr[i][j] << "," << yptr[i][j] << std::endl;
				}
				
			}
		}
	}

	template<typename T>
	void MatrixCopy(int copy_type,
					T*** X,
					T*** Y,
					int length,
					int start_xr = 0,
					int start_xc = 0,
					int start_yr = 0,
					int start_yc = 0) {
		auto& xptr = *(X);
		auto& yptr = *(Y);
		if (copy_type) { // from X to Y (X > Y)
			for (int i = start_xr; i < start_xr + length; i ++) {
				for (int j = start_xc; j < start_xc + length; j ++) {
					yptr[i - start_xr + start_yr][j - start_xc + start_yc] = xptr[i][j];
				}
			}
		} else {	// from Y to X (X > Y)
			for (int i = start_yr; i < start_yr + length; i ++) {
				for (int j = start_yc; j < start_yc + length; j ++) {
					xptr[i - start_yr + start_xr][j - start_yc + start_xc] = yptr[i][j];
				}
			}
		}
	}

	template<typename T>
	void CreateFoxMatrix(T*** A, int matrix_size, int init_type=1) {
		auto& aptr = *(A);
		aptr = new T* [matrix_size];
		for (int i = 0; i < matrix_size; i ++) {
			aptr[i] = new T [matrix_size];
			for (int j = 0; j < matrix_size; j ++) {
				if (init_type) {
					aptr[i][j] = rand() % 100;
				} else {
					aptr[i][j] = 0;
				}
			}
		}
	}

	template<typename T>
	void ReleaseFoxMatrix(T*** A, int matrix_size) {
		auto& aptr = *(A);
		for (int i = 0; i < matrix_size; i ++) {
			delete [] aptr[i];
		}
		delete [] aptr;
	}

	// 2. FOX Matrix Multiple

	void CreateFoxMatrixSingle(int& my_rank, int& psize, MPI_Comm my_comm, int matrix_size) {
		int** A = nullptr;
		int** B = nullptr;
		int** C = nullptr;
		srand(0);
		CreateFoxMatrix<int>(&A, matrix_size);
		CreateFoxMatrix<int>(&B, matrix_size);
		CreateFoxMatrix<int>(&C, matrix_size, 0);
		MatrixPrint<int>(&A, "A", matrix_size);
		MatrixPrint<int>(&B, "B", matrix_size);
		MatrixMult<int>(&A, &B, &C, matrix_size);
		MatrixPrint<int>(&C, "Normal-C", matrix_size);
		ReleaseFoxMatrix<int>(&A, matrix_size);
		ReleaseFoxMatrix<int>(&B, matrix_size);
		ReleaseFoxMatrix<int>(&C, matrix_size);
	}

	void CreateFoxMatrixMult(int& my_rank, int& psize, MPI_Comm my_comm) {
		int** A = nullptr;
		int** B = nullptr;
		int** C = nullptr;
		int** MA = nullptr;
		int** MB = nullptr;
		int sub_matrix_size = 256;
		int sqrt_size = (int)sqrt(psize);
		int matrix_size = sub_matrix_size * sqrt_size;
		srand(0);
		CreateFoxMatrix<int>(&A, matrix_size);
		CreateFoxMatrix<int>(&B, matrix_size);
		if (my_rank == 0) {
			MatrixPrint<int>(&A, "A", matrix_size);
			MatrixPrint<int>(&B, "B", matrix_size);
		}
		CreateFoxMatrix<int>(&C, matrix_size, 0);
		CreateFoxMatrix<int>(&MA, sub_matrix_size);
		CreateFoxMatrix<int>(&MB, sub_matrix_size);
		int row = (int)(my_rank / sqrt_size);
		int col = (int)(my_rank % sqrt_size);
		MPI_Comm row_comm;
		MPI_Comm col_comm;
		MPI_Comm_split(my_comm, row, col, &row_comm);
		MPI_Comm_split(my_comm, col, row, &col_comm);
		MatrixCopy<int>(1,
						&A,
						&MA,
						sub_matrix_size,
						row * sub_matrix_size,
						col * sub_matrix_size);
		MatrixCopy<int>(1,
						&B,
						&MB,
						sub_matrix_size,
						row * sub_matrix_size,
						col * sub_matrix_size);
		int** CopyMA = nullptr;
		int** MultC = nullptr;
		CreateFoxMatrix<int>(&CopyMA, sub_matrix_size);
		CreateFoxMatrix<int>(&MultC, sub_matrix_size);
		MPI_Barrier(my_comm);
		clock_t start_t = clock();
		// Fox multiple
		for (int i = 0; i < sqrt_size; i ++) {
			MatrixCopy<int>(1,
							&MA,
							&CopyMA,
							sub_matrix_size);
			for (int j = 0; j < sub_matrix_size; j ++) {
				MPI_Bcast(CopyMA[j], sub_matrix_size, MPI_INT, (i + row) % sqrt_size, row_comm);
			}
			MatrixMult<int>(&CopyMA, &MB, &MultC, sub_matrix_size);
			MatrixAdd<int>(&C, &MultC, &C, sub_matrix_size);
			// Shift MB
			MPI_Barrier(col_comm);
			for (int j = 0; j < sub_matrix_size; j ++) {
				MPI_Request request_s, request_r;
				MPI_Status status;
				int dest = (row == 0) ? sqrt_size - 1 : row - 1;
				int src = (row == (sqrt_size - 1)) ? 0 : row + 1;
				MPI_Isend(MB[j], sub_matrix_size, MPI_INT, dest, dest, col_comm, &request_s);
				MPI_Irecv(MB[j], sub_matrix_size, MPI_INT, src, row, col_comm, &request_r);
				MPI_Wait(&request_s, &status);
				MPI_Wait(&request_r, &status);
			}
		}
		// Collect matrix C
		if (my_rank != 0) {
			std::vector<MPI_Request> req_vec;
			for (int j = 0; j < sub_matrix_size; j ++) {
				MPI_Request request;
				// MPI_Status status;
				MPI_Isend(C[j], sub_matrix_size, MPI_INT, 0, my_rank, my_comm, &request);
				req_vec.push_back(request);
				// MPI_Wait(&request, &status);
			}
			MPI_Status status;
			for(auto iter = req_vec.begin(); iter != req_vec.end(); iter ++) {
				MPI_Wait(&(*iter), &status);
			}
		} else {
			std::vector<MPI_Request> req_vec;
			for (int i = 1; i < psize; i ++) {
				int row = (int)(i / sqrt_size);
				int col = (int)(i % sqrt_size);
				for (int j = 0; j < sub_matrix_size; j ++) {
					MPI_Request request;
					// MPI_Status status;
					MPI_Irecv(&C[row * sub_matrix_size + j][col * sub_matrix_size],
							 sub_matrix_size, MPI_INT, i, i, my_comm, &request);
					// MPI_Wait(&request, &status);
					req_vec.push_back(request);
				}	
			}
			MPI_Status status;
			for(auto iter = req_vec.begin(); iter != req_vec.end(); iter ++) {
				MPI_Wait(&(*iter), &status);
			}
		}
		clock_t end_t = clock();

		int use_t = (end_t - start_t);
		int* use_time_arr = new int[psize];
		MPI_Gather(&use_t, 1, MPI_INT, use_time_arr, 1, MPI_INT, 0, my_comm);

		if (my_rank == 0) {
			float mpi_sum_t = 0.0;
			for (int i = 0; i < psize; i ++) {
				mpi_sum_t += use_time_arr[i];
			}
			mpi_sum_t = mpi_sum_t / (float)psize;
			// MatrixPrint<int>(&C, "MPI-Fox-C", matrix_size);

			// Normal fox multpile
			int** NC = nullptr;
			CreateFoxMatrix<int>(&NC, matrix_size, 0);
			start_t = clock();
			MatrixMult<int>(&A, &B, &NC, matrix_size);
			end_t = clock();
			use_t = (end_t - start_t);
			MatrixCompare(&C, &NC, matrix_size);
			// MatrixPrint<int>(&NC, "Normal-Fox-C", matrix_size);
			std::cerr << matrix_size << "-" << sub_matrix_size << "-MPI-Fox-Avg-Time: " << mpi_sum_t << std::endl
					  << matrix_size << "-" << sub_matrix_size << "-Normal-Fox-Time: " << use_t << std::endl;
		}

		ReleaseFoxMatrix<int>(&A, matrix_size);
		ReleaseFoxMatrix<int>(&B, matrix_size);
		ReleaseFoxMatrix<int>(&C, matrix_size);
		ReleaseFoxMatrix<int>(&MA, sub_matrix_size);
		ReleaseFoxMatrix<int>(&MB, sub_matrix_size);
		ReleaseFoxMatrix<int>(&CopyMA, sub_matrix_size);
		ReleaseFoxMatrix<int>(&MultC, sub_matrix_size);
		delete [] use_time_arr;
	}

	// 3. Parameter Server
	void CreateParameterServer(int& my_rank, int& psize, MPI_Comm my_comm) {
		
		int pServerNum = 2;
		// int wServerNum = psize - pServerNum;

		int pColor = my_rank / pServerNum;
		int pKey = my_rank % pServerNum;

		MPI_Comm param_comm;
		MPI_Comm_split(my_comm, pColor, pKey, &param_comm);

		if (my_rank >= pServerNum) {
			srand(my_rank);
			while (true) {
				int sdata = rand() % 1000;
				int rdata = 0;
				int dest = my_rank % pServerNum;
				MPI_Status status;
				MPI_Request request;
#ifdef WINDOWS_PLATFORM
				int stime = rand() % 1000 + 2000;
				Sleep(stime);
#else
				int stime = rand() % 10;
				sleep(stime);
#endif // WINDOWS_PLATFROM
				std::cout << "* wServer " << my_rank << " Send to " << dest << " sData " << sdata << std::endl;
				MPI_Isend(&sdata, 1, MPI_INT, dest, my_rank, my_comm, &request);
				MPI_Wait(&request, &status);
				MPI_Irecv(&rdata, 1, MPI_INT, dest, my_rank, my_comm, &request);
				MPI_Wait(&request, &status);
				std::cout << "* wServer " << my_rank << " Recv from " << dest << " avg data " << rdata << std::endl;
			}
		} else {
			srand(my_rank);
			while (true) {
				int rNum = (psize - 1) / pServerNum - 1;
				if ((psize - 1) % pServerNum >= my_rank) {
					rNum++;
				}
				int* rData = new int[rNum];
				for (int i = 1; i <= rNum; i++) {
					int src = i * pServerNum + my_rank;
					MPI_Request request;
					MPI_Status status;
					MPI_Irecv(&rData[i], 1, MPI_INT, src, src, my_comm, &request);
					MPI_Wait(&request, &status);
					std::cout << "* pServer " << my_rank << " Recv from " << src << " rData " << rData[i] << std::endl;
				}
				for (int i = 2; i <= rNum; i++) {
					rData[1] += rData[i];
				}
				std::cout << "* pServer " << my_rank << " Compute data " << rData[1] << std::endl;
				int* gData = new int[pServerNum];
				MPI_Gather(&rData[1], 1, MPI_INT, gData, 1, MPI_INT, 0, param_comm);
				int avgData = 0;
				if (my_rank == 0) {
					std::cout << "* PServer rData: ";
					for (int i = 0; i < pServerNum; i ++) {
						std::cout << gData[i] << " ";
						avgData += gData[i];
					}
					avgData /= pServerNum;
					std::cout << std::endl;
				}
				MPI_Bcast(&avgData, 1, MPI_INT, 0, param_comm);
				std::cout << "* pServer " << my_rank << " Get avg data " << avgData << std::endl;
#ifdef WINDOWS_PLATFORM
				int stime = rand() % 1000 + 2000;
				Sleep(stime);
#else
				int stime = rand() % 10;
				sleep(stime);
#endif // WINDOWS_PLATFROM
				for (int i = 1; i <= rNum; i++) {
					int dest = i * pServerNum + my_rank;
					MPI_Request request;
					MPI_Status status;
					MPI_Isend(&avgData, 1, MPI_INT, dest, dest, my_comm, &request);
					MPI_Wait(&request, &status);
					std::cout << "* pServer " << my_rank << " Send to " << dest << " avg data " << avgData << std::endl;
				}
			}
		}
	}

	int NormalMod(int x, int y, int z) {
		// pow(x, y) % z
		int out = 1;
		for (int i = 0; i < y; i ++) {
			out *= x;
			out = out % z;
		}
		return out;
	}

	int FastMod(int x, int y, int z) {
		int ans = 1;
		int base = x % z;
		while (y) {
			if (y & 1) {
				ans = (ans * base) % z;
			}
			base = (base * base) % z;
			y >>= 1;
		}
		return ans;
	}

	// 4. MC Algorithm
	bool Btest(int& a, int& n, bool debug=false) {
		auto s = 0;
		auto t = n - 1;
		while (t % 2 == 0) {
			s ++;
			t = t / 2;
		}
		if (debug) {
			std::cout << "* n,a,t,s=" << n << "," << a << "," << t << "," << s << std::endl;
		}
		auto x = FastMod(a, t, n);
		if (debug) {
			std::cout << "* x=" << x << std::endl;
		}
		if (x == 1 || x == n - 1) {
			return true;
		}
		for (int i = 1; i < s; i ++) {
			x = FastMod(x, 2, n);
			if (debug) {
				std::cout << "* x=" << x << std::endl;
			}
			if (x == n - 1) {
				return true;
			}
		}
		return false;
	}

	bool MillerRabin(int& n, bool debug=false) {
		int a = rand() % ( n - 3) + 2;
		// std::cout << "* Rand choose: " << a << std::endl;
		return Btest(a, n, debug);
	}

	bool RepeatMillerRabin(int& n, int& repeat_num, bool debug=false) {
		// std::cout << "* MR: n,r=" << n << "," << repeat_num << std::endl;
		for (int i = 0; i < repeat_num; i ++) {
			if (!MillerRabin(n, debug)) {
				return false;
			}
		}
		return true;
	}

	bool PrimeJudge(int& n) {
		int sqrt_n = sqrt(n);
		// int sqrt_n = n;
		for (int i = 2; i <= sqrt_n; i ++) {
			if (n % i == 0) {
				return false;
			}
		}
		return true;
	}

	void CreateMonteCarloSingle(int& my_rank, int& psize, MPI_Comm my_comm) {
		int floor_n = 100000;
		std::vector<int> pj_vec;
		std::vector<int> mr_vec;

		srand(time(NULL));

		clock_t start_t = clock();
		for (int i = 5; i < floor_n; i += 2) {
			if (PrimeJudge(i)) {
				pj_vec.push_back(i);
			}
		}
		clock_t end_t = clock();

		// std::cout << "* PJ use time: " << (end_t - start_t) << "s" << std::endl;
		std::cerr << "1<PJ>Prime-size:" << pj_vec.size() << std::endl;
		std::cerr << "1<PJ>Avg-time-use(s):" << (end_t - start_t) << std::endl;
		start_t = clock();
		for (int i = 5; i < floor_n; i += 2) {
			int log_i = ceil(log(i));
			log_i = (int)ceil(sqrt(log_i));
			// int log_i = 10;
			if (RepeatMillerRabin(i, log_i)) {
				mr_vec.push_back(i);
			}
			// } else if (RepeatMillerRabin(i, log_i)) {
			// 	mr_vec.push_back(i);
			// }
		}
		end_t = clock();

		// std::cout << "* MR use time: " << (end_t - start_t) << "s" << std::endl;

		if (pj_vec.size() != mr_vec.size()) {
			std::cout << "* There is some difference" << std::endl;
		}
		// int j = 0;
		// int diff_cnt = 0;
		// for (int i = 0; i < pj_vec.size(); i ++) {
		// 	if (pj_vec[i] != mr_vec[j]) {
		// 		std::cout << "* Different " << pj_vec[i] << "," << mr_vec[j] << std::endl;
		// 		int log_i = ceil(log(pj_vec[i]));
		// 		RepeatMillerRabin(pj_vec[i], log_i, true);
		// 		diff_cnt ++;
		// 		// break;
		// 	} else {
		// 		j ++;
		// 	}
		// 	if (diff_cnt > 5) break;
		// }

		std::cerr << "1<MC>Prime-size:" << mr_vec.size() << std::endl;
		std::cerr << "1<MC>Avg-time-use(s):" << (end_t - start_t) << std::endl;
	}

	void CreateMonteCarloParallel(int& my_rank, int& psize, MPI_Comm my_comm) {

		typedef char bool_t;

		int floor_n = 100000000;
		int pn = ceil((float)(floor_n / 2) / (float)psize);
		int sum_pn = pn * psize;

		bool_t* prime_arr = new bool_t[pn];

		srand(time(NULL));

		clock_t start_t = clock();
		for (int i = 5 + my_rank * 2; i < floor_n; i += (psize * 2)) {
			int log_i = ceil(log(i));
			prime_arr[(i - 5) / (psize * 2)] = RepeatMillerRabin(i, log_i) ? '1' : '0';
		}
		
		bool_t* prime_sum_arr = new bool_t[sum_pn];
		MPI_Gather(prime_arr, pn, MPI_CHAR, prime_sum_arr, pn, MPI_CHAR, 0, my_comm);
		
		std::vector<int> prime_vec;
		if (my_rank == 0) {
			prime_vec.push_back(2);
			prime_vec.push_back(3);
			for (int p = 0; p < psize; p ++) {
				for (int i = 0; i < pn; i ++) {
					if (prime_sum_arr[p * pn + i] == '1') {
						prime_vec.push_back(5 + (p * 2) + (i * psize * 2));
					}
				}
			}
		}
		clock_t end_t = clock();

		int use_t = (int)(end_t - start_t);
		// std::cout << "* MR " << my_rank << " use time: " << use_t << "s" << std::endl;
		int* time_sum_arr = new int[psize];
		MPI_Gather(&use_t, 1, MPI_INT, time_sum_arr, 1, MPI_INT, 0, my_comm);
		if (my_rank == 0) {
			int time_sum = 0;
			for (int p = 0; p < psize; p ++) {
				time_sum += time_sum_arr[p];
			}
			// std::sort(prime_vec.begin(), prime_vec.end(), [] (int x, int y) {return x < y;});
			// std::cout << "* Prime array: ";
			// for (auto iter = prime_vec.begin(); iter != prime_vec.end(); iter ++) {
			// 	std::cout << *iter << ",";
			// }
			// std::cout << std::endl;
			std::cerr << psize << ">Prime-size:" << prime_vec.size() << std::endl;
			std::cerr << psize << ">Avg-time-use(s):" << (float)time_sum / (float)psize << std::endl;
		}

		delete [] prime_arr;
		delete [] prime_sum_arr;
	}

	enum init_t {TC, GS};

	// 15.1
	template<typename T>
	void CreateMatrix(T*** A, int row_size, int col_size, init_t init_type=init_t::TC) {
		auto& aptr = *(A);
		aptr = new T* [row_size];
		for (int i = 0; i < row_size; i ++) {
			aptr[i] = new T [col_size];
			for (int j = 0; j < col_size; j ++) {
				if (init_type == init_t::TC) {
					if (i != j) {
						int tmp = rand() % 10;
						aptr[i][j] = tmp / 8;
					} else {
						aptr[i][j] = 1;
					}
				} else if (init_type == init_t::GS) {
					aptr[i][j] = rand() % 20 + 1;
				}
			}
		}
	}

	template<typename T>
	void MatrixCopy(T*** A, T*** B, int row_size, int col_size) {
		auto& aptr = *(A);
		auto& bptr = *(B);
		for (int r = 0; r < row_size; r ++) {
			for (int c = 0; c < col_size; c ++) {
				bptr[r][c] = aptr[r][c];
			}
		}
	}

	template<typename T>
	void ReleaseMatrix(T*** A, int row_size) {
		auto& aptr = *(A);
		for (int r = 0; r < row_size; r ++) {
			delete [] aptr[r];
		}
		delete [] aptr;
	}

	template<typename T>
	void VectorPrint(T** X, std::string name, int vector_size) {
		auto& xptr = *(X);
		std::cout << "* VectorPrint " << name << std::endl;
		for (int i = 0; i < vector_size; i ++) {
			std::cout << xptr[i] << " ";
		}
		std::cout << std::endl;
	}

	template<typename T>
	void VectorCopy(T** X, T** Y, int vector_size) {
		auto& xptr = *(X);
		auto& yptr = *(Y);
		for (int i = 0; i < vector_size; i ++) {
			yptr[i] = xptr[i];
		}
	}

	template<typename T>
	bool VectorCompare(T** X, T** Y, int vector_size) {
		auto& xptr = *(X);
		auto& yptr = *(Y);
		for (int i = 0; i < vector_size; i ++) {
			auto diff = fabs(xptr[i] - yptr[i]);
			if (diff >= 0.003) {
				std::cout << "* Vector compare failed with " 
							<< xptr[i] << "," << yptr[i] << std::endl;
				return false;
			}
		}
		return true;
	}

	template<typename T>
	void CreateVector(T** X, int vector_size, init_t init_type=init_t::GS) {
		auto& xptr = *(X);
		xptr = new T[vector_size];
		for (int i = 0; i < vector_size; i ++) {
			if (init_type == init_t::GS) {
				xptr[i] = rand() % 20;
			} else {
				xptr[i] = 0;
			}
		}
	}

	template<typename T>
	void ReleaseVector(T** X) {
		auto& xptr = *(X);
		delete [] xptr;
	}

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

	void CreateTransitiveClosureParallel(int& my_rank, int& psize, MPI_Comm my_comm) {

		int** SA = nullptr;
		int** PA = nullptr;
		int** M = nullptr;
		int sub_matrix_size = 64;
		int matrix_size = sub_matrix_size * psize;
		srand(time(NULL));
		CreateMatrix<int>(&SA, matrix_size, matrix_size);
		CreateMatrix<int>(&PA, matrix_size, matrix_size);
		MatrixCopy<int>(&SA, &PA, matrix_size, matrix_size);
		CreateMatrix<int>(&M, matrix_size, matrix_size);
		int serial_use_t = 0;
		int parallel_use_t = 0;
		if (my_rank == 0) {
			// MatrixPrint<int>(&SA, "Initial-A", matrix_size);
			clock_t start_t = clock();
			RunTransitiveClosureSerial<int>(&SA, &M, matrix_size);
			clock_t end_t = clock();
			serial_use_t = (end_t - start_t);
			// MatrixPrint<int>(&SA, "Serial-A", matrix_size);
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

		clock_t start_t = clock();
		for (int i = 1; i <= log(matrix_size); i ++) {
			// 1. Root send sub_rows and sub_cols to each processor
			if (my_rank == 0) {
				std::vector<MPI_Request> req_vec;
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
				for (auto iter = req_vec.begin(); iter != req_vec.end(); iter ++) {
					MPI_Wait(&(*iter), &status);
				}
			} else {
				MPI_Request request;
				MPI_Status status;
				for (int r = 0; r < row_size; r ++) {
					MPI_Irecv(LocalA[r], matrix_size, MPI_INT, 0, my_rank, my_comm, &request);
					MPI_Wait(&request, &status);
				}
				for (int r = 0; r < matrix_size; r ++) {
					MPI_Irecv(LocalB[r], col_size, MPI_INT, 0, my_rank + psize, my_comm, &request);
					MPI_Wait(&request, &status);
				}
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
				for (int r = 0; r < matrix_size; r ++) {
					int src = (my_rank + 1) % psize;
					int dest = (my_rank == 0) ? psize - 1 : (my_rank - 1);
					MPI_Request request_s, request_r;
					// MPI_Status status;
					MPI_Isend(TmpB[r], col_size, MPI_INT, dest, my_rank, my_comm, &request_s);
					MPI_Irecv(LocalB[r], col_size, MPI_INT, src, src, my_comm, &request_r);
					req_vec_s.push_back(request_s);
					req_vec_r.push_back(request_r);
					// MPI_Wait(&request_s, &status);
					// MPI_Wait(&request_r, &status);
					// MPI_Status status;
					// MatrixCopy<int>(&LocalB, &TmpB, matrix_size, col_size);
					// if (my_rank == 0) {
					// 	MPI_Recv(LocalB[r], col_size, MPI_INT, src, src, my_comm, &status);
					// 	MPI_Send(TmpB[r], col_size, MPI_INT, dest, my_rank, my_comm);
					// } else {
					// 	MPI_Send(TmpB[r], col_size, MPI_INT, dest, my_rank, my_comm);
					// 	MPI_Recv(LocalB[r], col_size, MPI_INT, src, src, my_comm, &status);
					// }	
				}
				MPI_Status status;
				for (auto iter = req_vec_s.begin(); iter != req_vec_s.end(); iter ++) {
					MPI_Wait(&(*iter), &status);
				}
				for (auto iter = req_vec_r.begin(); iter != req_vec_r.end(); iter ++) {
					MPI_Wait(&(*iter), &status);
				}
			}
			// 3. Collect compute output
			if (my_rank != 0) {
				// Send sub_rows to root
				for (int r = 0; r < row_size; r ++) {
					MPI_Send(M[my_rank * row_size + r], matrix_size, MPI_INT, 0, my_rank, my_comm);
				}
			} else {
				// Recv sub_rows from other processors
				for (int p = 1; p < psize; p ++) {
					for (int r = 0; r < row_size; r ++) {
						MPI_Status status;
						MPI_Recv(M[p * row_size + r], matrix_size, MPI_INT, p, p, my_comm, &status);
					}	
				}
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
		clock_t end_t = clock();
		// Root show matrix A
		if (my_rank == 0) {
			// MatrixPrint<int>(&PA, "Parallel-A", matrix_size);
			MatrixCompare<int>(&SA, &PA, matrix_size);
		}
		// Compute use time
		parallel_use_t = (end_t - start_t);
		int* use_time_arr = new int[psize];
		MPI_Gather(&parallel_use_t, 1, MPI_INT, use_time_arr, 1, MPI_INT, 0, my_comm);
		if (my_rank == 0) {
			float mpi_sum_t = 0.0;
			for (int i = 0; i < psize; i ++) {
				mpi_sum_t += use_time_arr[i];
			}
			mpi_sum_t = mpi_sum_t / (float)psize;
			std::cerr << psize << "\t" << sub_matrix_size << "\t" << serial_use_t << "\t" << mpi_sum_t << std::endl;
		}

		ReleaseMatrix<int>(&SA, matrix_size);
		ReleaseMatrix<int>(&PA, matrix_size);
		ReleaseMatrix<int>(&M, matrix_size);
		ReleaseMatrix<int>(&LocalA, row_size);
		ReleaseMatrix<int>(&LocalB, matrix_size);
		ReleaseMatrix<int>(&TmpB, matrix_size);
	}

	// 19.1 Gauss Elimination
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
			// MatrixPrint<T>(&aptr, "Swith-A", matrix_size);
			// VectorPrint<T>(&bptr, "Swith-B", matrix_size);
			// Div current line
			for (int c = k + 1; c < col_size; c ++) {
				aptr[k][c] = aptr[k][c] / aptr[k][k];
			}
			bptr[k] = bptr[k] / aptr[k][k];
			aptr[k][k] = 1;
			// MatrixPrint<T>(&aptr, "DIV-A", matrix_size);
			// VectorPrint<T>(&bptr, "DIV-B", matrix_size);
			// Start eliminate
			for (int r = k + 1; r < row_size; r ++) {
				for (int c = k + 1; c < col_size; c ++) {
					aptr[r][c] -= aptr[r][k] * aptr[k][c];
				}
				bptr[r] -= aptr[r][k] * bptr[k];
				aptr[r][k] = 0;
			}
			// MatrixPrint<T>(&aptr, "GSE-A", row_size, col_size);
			// VectorPrint<T>(&bptr, "GSE-B", row_size);
		}
		// MatrixPrint<T>(&aptr, "GSE-A", row_size, col_size);
		// VectorPrint<T>(&bptr, "GSE-B", row_size);
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
		// VectorPrint<T>(&xptr, "Serial-X", col_size);
		ReleaseVector<T>(&SortX);
	}

	void CreateGaussEliminParallel(int& my_rank, int& psize, MPI_Comm my_comm) {
		float** SA = nullptr;
		float** PA = nullptr;
		float* SB = nullptr;
		float* PB = nullptr;
		float* SX = nullptr;
		float* PX = nullptr;
		int sub_matrix_size = 64;
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
			// MatrixPrint<float>(&SA, "Initial-A", row_size, col_size);
			// VectorPrint<float>(&SB, "Initial-B", row_size);
			RunGaussEliminSerial<float>(&SA, &SB, &SX, row_size, col_size);
			ReleaseMatrix<float>(&SA, row_size);
			ReleaseVector<float>(&SB);
			// ReleaseVector<float>(&SX);
		}
		MPI_Barrier(my_comm);

		float** LocalA = nullptr;
		float* LocalB = nullptr;
		float* GlobalF = nullptr;
		CreateMatrix<float>(&LocalA, sub_matrix_size, col_size);
		CreateVector<float>(&LocalB, sub_matrix_size);
		CreateVector<float>(&GlobalF, col_size + 1);

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
				MPI_Alltoall(lmax_all_put, 4, MPI_FLOAT, lmax_all_get, 4, MPI_FLOAT, my_comm);
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
						MPI_Recv(getTmpA, col_size + 1, MPI_FLOAT, (int)lmax[3], (int)lmax[3], my_comm, &status);
						MPI_Send(putTmpA, col_size + 1, MPI_FLOAT, (int)lmax[3], my_rank, my_comm);
						VectorCopy<float>(&getTmpA, &LocalA[s], col_size);
						LocalB[s] = getTmpA[col_size];
					}
					if (my_rank == (int)lmax[3]) {
						// send
						VectorCopy<float>(&LocalA[(int)lmax[1]], &putTmpA, col_size);
						putTmpA[col_size] = LocalB[(int)lmax[1]];
						MPI_Send(putTmpA, col_size + 1, MPI_FLOAT, p, my_rank, my_comm);
						MPI_Recv(getTmpA, col_size + 1, MPI_FLOAT, p, p, my_comm, &status);
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
				MPI_Bcast(GlobalF, col_size + 1, MPI_FLOAT, p, my_comm);
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
					MPI_Bcast(&PX[s * psize + p], 1, MPI_FLOAT, p, my_comm);
					for (int ss = 0; ss < sub_matrix_size; ss ++) {
						sum[ss] += LocalA[ss][s * psize + p] * PX[s * psize + p];
					}
				} else {
					MPI_Bcast(&PX[s * psize + p], 1, MPI_FLOAT, p, my_comm);
					int end_row = (my_rank > p) ? s - 1 : s;
					for (int ss = 0; ss <= end_row; ss ++) {
						sum[ss] += LocalA[ss][s * psize + p] * PX[s * psize + p];
					}
					// std::cout << "* PX=" << PX[s * psize + p] << " end_row:" << end_row 
					// 	<< " sum=" << sum[end_row - 1] << std::endl;
				}
				MPI_Barrier(my_comm);
			}
		}
		// VectorPrint<float>(&PX, "Parallel-X", col_size);
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
			// VectorPrint<float>(&PX, "Parallel-X", col_size);
			if (VectorCompare<float>(&SX, &PX, col_size)) {
				std::cout << "* Serial vs Parallel pass !" << std::endl;
			}
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

	// 22.1 FFT
	typedef struct {
		double r;
		double i;
	} complex_t;

	complex_t ComplexMult(complex_t comp_a, complex_t comp_b) {
		complex_t comp_c;
		comp_c.r = comp_a.r * comp_b.r - comp_a.i * comp_b.i;
		comp_c.i = comp_a.r * comp_b.i + comp_a.i * comp_b.r;
		return comp_c;
	}

	complex_t ComplexAdd(complex_t comp_a, complex_t comp_b) {
		complex_t comp_c;
		comp_c.r = comp_a.r + comp_b.r;
		comp_c.i = comp_a.i + comp_b.i;
		return comp_c;
	}

	complex_t ComplexSub(complex_t comp_a, complex_t comp_b) {
		complex_t comp_c;
		comp_c.r = comp_a.r - comp_b.r;
		comp_c.i = comp_a.i - comp_b.i;
		return comp_c;
	}

	void ComplexInit(complex_t** comp_a, int vec_size, int init_type=0) {
		auto& vec_a = *(comp_a);
		for (int c = 0; c < vec_size; c ++) {
			if (init_type == 1) {
				vec_a[c].r = rand() % 10 + 1;
				vec_a[c].i = 0;
			} else {
				vec_a[c].r = cos((c * 2 * PI) / vec_size);
				vec_a[c].i = sin((c * 2 * PI) / vec_size);
			}
		}
	}

	void ComplexCopy(complex_t** comp_a, complex_t** comp_b, int vec_size) {
		auto& vec_a = *(comp_a);
		auto& vec_b = *(comp_b);
		for (int c = 0; c < vec_size; c ++) {
			vec_b[c].r = vec_a[c].r;
			vec_b[c].i = vec_a[c].i;
		}
	}

	void ComplexPrint(complex_t** comp_a, std::string name, int vec_size) {
		auto& vec_a = *(comp_a);
		std::cout << "* ComplexPrint " << name << std::endl;
		for (int c = 0; c < vec_size; c ++) {
			std::cout << "(" << vec_a[c].r << "," << vec_a[c].i << ") ";
		}
		std::cout << std::endl;
	}

	bool ComplexCompare(complex_t** comp_a, complex_t** comp_b, int vec_size) {
		auto& vec_a = *(comp_a);
		auto& vec_b = *(comp_b);
		for (int c = 0; c < vec_size; c ++) {
			if ((vec_a[c].r != vec_b[c].r) || (vec_a[c].i != vec_b[c].i)) {
				std::cout << "* Vector compare failed with " 
					<< "(" << vec_a[c].r << "," << vec_a[c].i << ") "
					<< "(" << vec_b[c].r << "," << vec_b[c].i << ") " << std::endl;
				return false;
			}
		}
		return true;
	}

	void CreateFftParallel(int& my_rank, int& psize, MPI_Comm my_comm) {

		int n = 1024;

		complex_t* comp_vec_a = new complex_t[n];
		complex_t* comp_vec_pa = new complex_t[n];
		complex_t* comp_vec_b = new complex_t[n];
		complex_t* comp_vec_w = new complex_t[n * 2];
		ComplexInit(&comp_vec_a, n, 1);
		ComplexCopy(&comp_vec_a, &comp_vec_pa, n);
		ComplexInit(&comp_vec_w, 2 * n);

		if (my_rank == 0) {
			// ComplexPrint(&comp_vec_a, "Init-A", n);
			// ComplexPrint(&comp_vec_w, "Init-W", 2 * n);
			for (int s = log(n) / log(2) - 1; s >= 0; s --) {
				int part = pow(2, s);
				int heap = n / part;
				int index_w = heap / 2;
				// std::cout << "* Part " << s << " " << part << std::endl;
				for (int c = 0; c < n; c ++) {
					int mod_p = c % part;
					int mod_2p = c % (2 * part);
					if (mod_p == mod_2p) {
						int sub_index_w = c % part;
						complex_t tmp_a = comp_vec_a[c];
						comp_vec_a[c] = ComplexAdd(comp_vec_a[c], comp_vec_a[c + part]);
						comp_vec_a[c + part] = ComplexSub(tmp_a, comp_vec_a[c + part]);
						comp_vec_a[c + part] = ComplexMult(comp_vec_a[c + part], comp_vec_w[sub_index_w * index_w]);
					}
				}
			}
			// ComplexPrint(&comp_vec_a, "Serial-A", n);			
		}
		MPI_Barrier(my_comm);
		int m = n / psize;
		// 1.
		for (int s = log(n) / log(2) - 1; s >= log(m) / log(2); s --) {
			int part = pow(2, s);
			int heap = n / part;
			int index_w = heap / 2;
			int startPos = my_rank * m;
			int mod_p = startPos % part;
			int mod_2p = startPos % (2 * part);
			// std::cout << "* Rank " << my_rank << " Part " << s << " " << part << std::endl;
			if (mod_p == mod_2p) {
				MPI_Status status;
				int backPos = (my_rank * m + part);
				MPI_Send(&comp_vec_pa[startPos], m * 2, MPI_DOUBLE, backPos / m, my_rank, my_comm);
				MPI_Recv(&comp_vec_pa[backPos], m * 2, MPI_DOUBLE, backPos / m, backPos / m, my_comm, &status);
				for (int c = startPos; c < startPos + m; c ++) {
					comp_vec_pa[c] = ComplexAdd(comp_vec_pa[c], comp_vec_pa[c + part]);					
				}
				
			} else {
				MPI_Status status;
				int frontPos = (my_rank * m - part);
				MPI_Recv(&comp_vec_pa[frontPos], m * 2, MPI_DOUBLE, frontPos / m, frontPos / m, my_comm, &status);
				MPI_Send(&comp_vec_pa[startPos], m * 2, MPI_DOUBLE, frontPos / m, my_rank, my_comm);
				for (int c = startPos; c < startPos + m; c ++) {
					int sub_index_w = c % part;
					comp_vec_pa[c] = ComplexSub(comp_vec_pa[c - part], comp_vec_pa[c]);
					comp_vec_pa[c] = ComplexMult(comp_vec_pa[c], comp_vec_w[sub_index_w * index_w]);					
				}
			}
		}
		// 2.
		for (int s = log(m) / log(2) - 1; s >= 0; s --) {
			int part = pow(2, s);
			int heap = n / part;
			int index_w = heap / 2;
			int startPos = my_rank * m;
			// std::cout << "* Rank " << my_rank << " Part " << s << " " << part << std::endl;
			for (int c = startPos; c < startPos + m; c ++) {
				int mod_p = c % part;
				int mod_2p = c % (2 * part);
				if (mod_p == mod_2p) {
					int sub_index_w = c % part;
					complex_t tmp_a = comp_vec_pa[c];
					comp_vec_pa[c] = ComplexAdd(comp_vec_pa[c], comp_vec_pa[c + part]);
					comp_vec_pa[c + part] = ComplexSub(tmp_a, comp_vec_pa[c + part]);
					comp_vec_pa[c + part] = ComplexMult(comp_vec_pa[c + part], comp_vec_w[sub_index_w * index_w]);
				}
			}
		}
		// 3.
		MPI_Gather(&comp_vec_pa[my_rank * m], m * 2, MPI_DOUBLE, comp_vec_b, m * 2, MPI_DOUBLE, 0, my_comm);
		if (my_rank == 0) {
			// ComplexPrint(&comp_vec_b, "Parallel-A", n);
			if (ComplexCompare(&comp_vec_a, &comp_vec_b, n)) {
				std::cout << "* Serial vs. Parallel pass !" << std::endl;
			}
		}

		delete [] comp_vec_a;
		delete [] comp_vec_pa;
		delete [] comp_vec_b;
		delete [] comp_vec_w;
	}
	
}

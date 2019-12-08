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
		MPI_Status status_r, status_s;
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

		// Step
		for (int s = 0; pow(2, s + 1) <= n; s++) {
			int div = pow(2, s + 1);
			int sub = pow(2, s);
			
			for (int i = 0; i < n; i++) {
				if (my_rank == i && (i % sub) == 0) {
					if ((i % div) != 0) {
						int dest = my_rank - sub;
						std::cout << "* Rank " << my_rank << " Send to " << dest << std::endl;
						MPI_Send(&sumArray[my_rank], 1, MPI_INT, dest, s, my_comm);
					} else {
						int src = my_rank + sub;
						MPI_Status status;
						MPI_Recv(&sumArray[src], 1, MPI_INT, src, s, my_comm, &status);
						sumArray[my_rank] += sumArray[src];
						std::cout << "* Rank " << my_rank << " Recv from " << src << std::endl;
					}
				}
			}
			MPI_Barrier(my_comm);
		}

		if (my_rank == 0) {
			if (sumArray[0] == check_sum) {
				std::cout << "* Rank: " << my_rank << " Sum: " << sumArray[my_rank] << std::endl;
			}
			else {
				std::cerr << "* Rank: " << my_rank << " Sum Check failed !" << std::endl;
			}
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
	void MatrixPrint(T*** X, std::string name, int matrix_size) {
		auto& xptr = *(X);
		std::cout << "* MatrixPrint " << name << std::endl;
		for (int i = 0; i < matrix_size; i ++) {
			for (int j = 0; j < matrix_size; j ++) {
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
		// MPI_Barrier(my_comm);
		// for (int i = 0; i < psize; i ++) {
		// 	MPI_Barrier(my_comm);
		// 	if (i == my_rank) {
		// 		std::string matrix_str("OutC_" + std::to_string(my_rank));
		// 		MatrixPrint<int>(&C, matrix_str, sub_matrix_size);
		// 	}
			
		// }
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
		int wServerNum = psize - pServerNum;

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
	bool Btest(int& a, int& n) {
		auto s = 0;
		auto t = n - 1;
		while (t % 2 == 0) {
			s ++;
			t = t / 2;
		}
		// std::cout << "* n,a,t,s=" << n << "," << a << "," << t << "," << s << std::endl;
		auto x = FastMod(a, t, n);
		// std::cout << "* x=" << x << std::endl;
		if (x == 1 || x == n - 1) {
			return true;
		}
		for (int i = 1; i < s; i ++) {
			x = FastMod(x, 2, n);
			// std::cout << "* x=" << x << std::endl;
			if (x == n - 1) {
				return true;
			}
		}
		return false;
	}

	bool MillerRabin(int& n) {
		int a = rand() % ( n - 3) + 2;
		// std::cout << "* Rand choose: " << a << std::endl;
		return Btest(a, n);
	}

	bool RepeatMillerRabin(int& n, int& repeat_num) {
		// std::cout << "* MR: n,r=" << n << "," << repeat_num << std::endl;
		for (int i = 0; i < repeat_num; i ++) {
			if (!MillerRabin(n)) {
				return false;
			}
		}
		return true;
	}

	bool PrimeJudge(int& n) {
		// int sqrt_n = sqrt(n);
		int sqrt_n = n;
		for (int i = 2; i < sqrt_n; i ++) {
			if (n % i == 0) {
				return false;
			}
		}
		return true;
	}

	void CreateMonteCarloSingle(int& my_rank, int& psize, MPI_Comm my_comm) {
		int floor_n = 100000000;
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
			// int log_i = 10;
			if (RepeatMillerRabin(i, log_i)) {
				mr_vec.push_back(i);
			}
		}
		end_t = clock();

		// std::cout << "* MR use time: " << (end_t - start_t) << "s" << std::endl;

		if (pj_vec.size() != mr_vec.size()) {
			std::cout << "* There is some difference" << std::endl;
		}

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
}

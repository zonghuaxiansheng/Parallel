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

// #define NUM_THREADS 4

namespace ustc_parallel {

double sstart_t = 0.0, send_t = 0.0, pstart_t = 0.0, pend_t = 0.0, trans_t = 0.0;

void CreatePrimeParallel(int& n, int& my_rank, int& psize, MPI_Comm my_comm) {

    typedef char bool_t;
    std::vector<int> prime_svec;
    std::vector<int> prime_pvec;
    srand(time(NULL));

    if (my_rank == 0) {
        sstart_t = MPI_Wtime();
        // prime_svec.push_back(2);
        // prime_svec.push_back(3);
        // for (int i = 5; i < n; i += 2) {
        //     int sqrt_n = sqrt(i);
        //     int is_prime = true;
        //     for (int j = 2; j <= sqrt_n; j ++) {
        //         if (i % j == 0) {
        //             is_prime = false;
        //             break;
        //         }
        //     }
        //     if (is_prime) {
        //         prime_svec.push_back(i);
        //     }
        // }
        send_t = MPI_Wtime();
    }
    MPI_Barrier(my_comm);
    int m = n / (psize * 2);
    int num_threads = 1;
    bool_t* prime_b = new bool_t[m];
    pstart_t = MPI_Wtime();
    // OMP
    #pragma omp parallel shared(num_threads)
    {
        num_threads = omp_get_num_threads();
        int my_id = omp_get_thread_num();
        int sub_m = m / num_threads;
        #pragma omp parallel for
        for (int i = 0; i < sub_m ; i ++) {
            int k = 5 + (my_rank * 2) + (i * psize * 2) + (sub_m * my_id * 2);
            int sqrt_n = sqrt(k);
            bool is_prime = true;
            for (int j = 2; j <= sqrt_n; j ++) {
                if (k % j == 0) {
                    is_prime = false;
                    break;
                }
            }
            prime_b[i + sub_m * my_id] = is_prime ? '1' : '0';     
        }
    }

    bool_t* prime_arr = new bool_t[n / 2];
    double tmp_t = MPI_Wtime();
    MPI_Gather(prime_b, m, MPI_CHAR, prime_arr, m, MPI_CHAR, 0, my_comm);
    trans_t += (MPI_Wtime() - tmp_t);
    pend_t = MPI_Wtime();
    if (my_rank == 0) {
        prime_pvec.push_back(2);
        prime_pvec.push_back(3);
        for (int p = 0; p < psize; p ++) {
            for (int i = 0; i < m; i ++) {
                int k = 5 + (p * 2) + (i * psize * 2);
                if (prime_arr[p * m + i] == '1' && k < n) {
                    prime_pvec.push_back(k);
                }
            }
        }
        if (prime_svec.size() == prime_pvec.size()) {
            std::cout << "* Serial vs. Parallel pass ! " << prime_svec.size() << "," << prime_pvec.size() << std::endl;
        } else {
            std::cout << "* Serial vs. Parallel failed ! " << prime_svec.size() << "," << prime_pvec.size() << std::endl;
        }
        std::cerr << n << "\t" << psize << "\t" << num_threads << "\t" << (send_t - sstart_t) << "\t" 
            << (pend_t - pstart_t) << "\t" << trans_t << "\t" << (pend_t - pstart_t - trans_t) << std::endl;
    }

    delete [] prime_b;
    delete [] prime_arr;
}

}

int main(int argc, char* argv[]) {

    int my_rank, psize;
    int n = 100000000;
    
    if (argc == 2) {
        n = atoi(argv[1]);
    }

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &psize);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);   

    ustc_parallel::CreatePrimeParallel(n, my_rank, psize, MPI_COMM_WORLD);

    MPI_Finalize();

	return 0;
}
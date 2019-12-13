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

#define PI    3.1415926535897932

namespace ustc_parallel {

double sstart_t = 0.0, send_t = 0.0, pstart_t = 0.0, pend_t = 0.0, trans_t = 0.0;

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

void CreateFftParallel(int& n, int& my_rank, int& psize, MPI_Comm my_comm) {

    complex_t* comp_vec_a = new complex_t[n];
    complex_t* comp_vec_pa = new complex_t[n];
    complex_t* comp_vec_b = new complex_t[n];
    complex_t* comp_vec_w = new complex_t[n * 2];
    ComplexInit(&comp_vec_a, n, 1);
    ComplexCopy(&comp_vec_a, &comp_vec_pa, n);
    ComplexInit(&comp_vec_w, 2 * n);

    if (my_rank == 0) {
        sstart_t = MPI_Wtime();
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
        send_t = MPI_Wtime();
        // ComplexPrint(&comp_vec_a, "Serial-A", n);			
    }
    MPI_Barrier(my_comm);
    pstart_t = MPI_Wtime();
    int m = n / psize;
    // 1.
    for (int s = log(n) / log(2) - 1; s >= log(m) / log(2); s --) {
        int part = pow(2, s);
        int heap = n / part;
        int index_w = heap / 2;
        int startPos = my_rank * m;
        int mod_p = startPos % part;
        int mod_2p = startPos % (2 * part);
        if (mod_p == mod_2p) {
            MPI_Status status;
            int backPos = (my_rank * m + part);
            double tmp_t = MPI_Wtime();
            MPI_Send(&comp_vec_pa[startPos], m * 2, MPI_DOUBLE, backPos / m, my_rank, my_comm);
            MPI_Recv(&comp_vec_pa[backPos], m * 2, MPI_DOUBLE, backPos / m, backPos / m, my_comm, &status);
            trans_t += (MPI_Wtime() - tmp_t);
            for (int c = startPos; c < startPos + m; c ++) {
                comp_vec_pa[c] = ComplexAdd(comp_vec_pa[c], comp_vec_pa[c + part]);					
            }
            
        } else {
            MPI_Status status;
            int frontPos = (my_rank * m - part);
            double tmp_t = MPI_Wtime();
            MPI_Recv(&comp_vec_pa[frontPos], m * 2, MPI_DOUBLE, frontPos / m, frontPos / m, my_comm, &status);
            MPI_Send(&comp_vec_pa[startPos], m * 2, MPI_DOUBLE, frontPos / m, my_rank, my_comm);
            trans_t += (MPI_Wtime() - tmp_t);
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
    pend_t = MPI_Wtime();
    MPI_Gather(&comp_vec_pa[my_rank * m], m * 2, MPI_DOUBLE, comp_vec_b, m * 2, MPI_DOUBLE, 0, my_comm);
    if (my_rank == 0) {
        // ComplexPrint(&comp_vec_b, "Parallel-A", n);
        if (ComplexCompare(&comp_vec_a, &comp_vec_b, n)) {
            std::cout << "* Serial vs. Parallel pass !" << std::endl;
        }
        std::cerr << n << "\t" << psize << "\t" << (send_t - sstart_t) << "\t" 
            << (pend_t - pstart_t) << "\t" << trans_t << "\t" << (pend_t - pstart_t - trans_t) << std::endl;
    }

    delete [] comp_vec_a;
    delete [] comp_vec_pa;
    delete [] comp_vec_b;
    delete [] comp_vec_w;
}
}

int main(int argc, char* argv[]) {

    int my_rank, psize;
    int n = 1024;
    
    if (argc == 2) {
        n = atoi(argv[1]);
    }

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &psize);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);   

    ustc_parallel::CreateFftParallel(n, my_rank, psize, MPI_COMM_WORLD);

    MPI_Finalize();

	return 0;
}
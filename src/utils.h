#ifndef _UTILS_H_
#define _UTILS_H_
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

namespace ustc_parallel {

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

enum init_t {TC, GS};

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
}
#endif

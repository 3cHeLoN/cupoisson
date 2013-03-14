/* 
 *  utils.cu - this file is part of CuPoisson
 *
 *  Copyright © 2011-2013, Folkert Bleichrodt 
 *	
 *  CuPoisson is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  CuPoisson is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with CuPoisson.  If not, see <http://www.gnu.org/licenses/>.
 */


/*
 * These functions are not necessary for running the code. However, they
 * can be used to simplify error checking and safeguard CUDA kernels.
 */

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cublas.h>
#include <string.h>
#include "utils.cuh"

#define BLOCK_SIZE 16
#define THREADSPB (BLOCK_SIZE * BLOCK_SIZE)
#define IDX(i,j,ld) ((i-1)*ld+(j-1))


/* initialize a vector to a constant value */
__global__ void cuFillArray(unsigned int n, real *dest, const real value)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= n) return;

    dest[idx] = value;
}

/* error checking */
void checkCudaError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if(cudaSuccess != err) {
    fprintf(stderr, "Cuda error: %s: %s.\n", msg, 
                        cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }                         
}

void writeFile(unsigned int gridSizeX, unsigned int gridSizeY, const real *x, const char *filename)
{
    unsigned int n = (gridSizeX-2) * (gridSizeY-2);
    real *temp = (real *)malloc(n*sizeof(real));

    if (cublasGetVector(n, sizeof(real), x, 1, temp, 1) != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "Error: cublas failed!\n");
        exit(EXIT_FAILURE);
    }
    
    FILE *fp = fopen(filename, "w");

    // leading dimensions gridsize
    int ld = gridSizeY - 2;
    
    int i, j;
    for (j = 0; j < gridSizeY; j++) {
        for (i = 0; i < gridSizeX; i++) {
            if (i == 0 || j == 0 || i == gridSizeX - 1 || j == gridSizeY - 1)
                fprintf(fp, "0 ");
            else
                fprintf(fp, "%f ", temp[IDX(i,j,ld)]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    free(temp);
}

void writeBinaryFile(unsigned int gridSizeX, unsigned int gridSizeY, const real *x, const char *filename)
{
    size_t FSIZE = sizeof(real);

    unsigned int n = (gridSizeX-2) * (gridSizeY-2);
    real *temp = (real*)malloc(n*FSIZE);

    if (cublasGetVector(n, FSIZE, x, 1, temp, 1) != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "Error: cublas failed!\n");
        exit(EXIT_FAILURE);
    }

    FILE *fp = fopen(filename, "wb");

    // header
    int dimensions[2] = {gridSizeX, gridSizeY};

    fwrite(dimensions, sizeof(int), 2, fp);

    real *zeros;

    // boundary conditions
    zeros = (real*)malloc(gridSizeY*FSIZE);
    int i;

    for (i = 0; i < gridSizeY; i++) {
        zeros[i] = 0.0;
    }

    fwrite(zeros, FSIZE, gridSizeY, fp);

    for (i = 0; i < gridSizeX-2; i++) {
        // boundary condition
        fwrite(zeros, FSIZE, 1, fp);
        // one column
        fwrite(temp+i*(gridSizeY-2),  FSIZE, gridSizeY-2, fp);
        // boundary condition
        fwrite(zeros, FSIZE, 1, fp);
    }

    fwrite(zeros, FSIZE, gridSizeY, fp);

    // clean up
    fclose(fp);
    free(temp);
    free(zeros);
}


/**
 * This function is a wrapper for the CUDA kernel
 * cuFillArray, to use in C code without loading CUDA specific
 * libraries.
 * 
 * @param dest,  device pointer for destination
 * @param value, the value to set for each element
 * @param count, number of elements to fill
 */
void fillArray(real *dest, const real value, unsigned int count)
{
    // determine number of blocks needed
    int nBlocks  = (count + THREADSPB-1)/THREADSPB;

    // call CUDA kernel
    cuFillArray<<<nBlocks, THREADSPB>>>(count, dest, value);
}


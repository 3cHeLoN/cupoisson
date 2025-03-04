/* 
 *  poisson.cu - this file is part of CuPoisson
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
 *  Poisson solver using CUFFT
 *  -------------------------
 * 
 *  Solves:   __
 *           -\/ u(x,y) = z(x,y);
 *  with
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cublas.h>
#include <cufft.h>
#include "poisson.cuh"
#include "precision.h"

// thread-block dimensions
#define BLOCK_SIZE 16
#define THREADSPB (BLOCK_SIZE * BLOCK_SIZE)

#define PI 3.1415926535897932384626433832795028841971693993751

/**
* @param n, size of input data (grid has size n x n)
* @param u, vector of unknowns
* @param z, right-hand side
* @return exit code for error checking
*/
int cuPoisson(int n, real *u, real *z)
{
    int m = 2*(n+1);
    
    // thread block dimensions
    dim3 dimBlock, dimGrid;
    dimBlock.x = dimBlock.y = BLOCK_SIZE;
    dimGrid.x = dimGrid.y = (n + BLOCK_SIZE-1)/BLOCK_SIZE;
    int nBlocks = (n + THREADSPB-1)/THREADSPB;
    
    // allocate a buffer
    // for preparing DST-I using FFTs
    real *buffer;
    cudaMalloc((void **) &buffer, n*m*sizeof(real));

    // setup CUFFT
    cufftHandle plan;

    // realFFT, half the size: m/2
    #ifdef DOUBLE_PRECISION
        cufftPlan1d(&plan, n+1, CUFFT_Z2Z, n);
    #else
        cufftPlan1d(&plan, n+1, CUFFT_C2C, n);
    #endif
    cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_NATIVE);



    //*DST-I(A) row-wise
    init_columns<<<nBlocks, THREADSPB>>>(n, m, buffer);

    // fill buffer to prepare batch FFTs
    copy_flip<<<dimGrid, dimBlock>>>(n, m, z, buffer);
    // cast real* into complex* of half the size
    cufftExec(plan, (complex *)buffer, (complex *)buffer, CUFFT_FORWARD);

    // extract realFFT and DST-I, also transpose
    extract_rfft<<<dimGrid, dimBlock>>>(n, m, (complex *)buffer, u);


    //*DST-I(B) column-wise (by transposing rows)
    init_columns<<<nBlocks, THREADSPB>>>(n, m, buffer);
    copy_flip<<<dimGrid, dimBlock>>>(n, m, u, buffer);
    cufftExec(plan, (complex *)buffer, (complex *)buffer, CUFFT_FORWARD);

    // again extract realFFT and transpose
    extract_rfft<<<dimGrid, dimBlock>>>(n, m, (complex *)buffer, u);

    // scale by eigenvalue and prepare for 2nd DST
    init_columns<<<nBlocks, THREADSPB>>>(n, m, buffer);
    scale_copy_flip<<<dimGrid, dimBlock>>>(n, m, u, buffer);



    //*DST-I(A) rows
    cufftExec(plan, (complex *)buffer, (complex *)buffer, CUFFT_FORWARD);
    extract_rfft<<<dimGrid, dimBlock>>>(n, m, (complex *)buffer, u);

    //*DST-I(B) columns
    init_columns<<<nBlocks, THREADSPB>>>(n, m, buffer);
    copy_flip<<<dimGrid, dimBlock>>>(n, m, u, buffer);
    cufftExec(plan, (complex *)buffer, (complex *)buffer, CUFFT_FORWARD);

    // transpose one more time and apply scaling
    // of the DST-I, copy back in place
    extract_scale<<<dimGrid, dimBlock>>>(n, m, 1.0/((n+1.0)*(n+1.0)), (complex *)buffer, u);

    // clean up CUFFT
    cufftDestroy(plan);
    cudaFree(buffer);

    return 0;
}

/**
 * Wrapper function for CUFFT
 */
cufftResult cufftExec(cufftHandle plan, complex *idata, complex *odata, int direction)
{
    cufftResult result;
    #ifdef DOUBLE_PRECISION
        result = cufftExecZ2Z(plan, idata, odata, direction);
    #else
        result = cufftExecC2C(plan, idata, odata, direction);
    #endif
    
    return result;
}

/*
* This function extracts the real FFT from the output of the complex FFT, but takes
* only the imaginary part, which is needed for the DST-I. The DST-I is reconstructed
* and the matrix is transposed to prepare for the columnwise DST-I.
*/
__global__ void extract_rfft(int n, int m, complex *src, real *dst)
{
    __shared__ real data[BLOCK_SIZE][BLOCK_SIZE+1];
    
    // column number
    int x = blockIdx.x * blockDim.x + threadIdx.x + 1;
    // row number
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    
    // global index (for complex type)
    int idx_l = y*(n+1)+x;
    int idx_r = (y+1)*(n+1)-x;

    // read part of data to shared mem
    if (x < n+1 && y < n)
        data[threadIdx.y][threadIdx.x] = src[idx_l].y;

    __syncthreads();

    // extract RFFT to shared mem
    if (x < n+1 && y < n) {
        data[threadIdx.y][threadIdx.x] -= 0.5*((data[threadIdx.y][threadIdx.x]+src[idx_r].y)*(1.0 + sin((2.0*PI*x)/m)) + (src[idx_l].x-src[idx_r].x)*cos((2.0*PI*x)/m));
    }
    __syncthreads();

    // transposing matrix
    x = blockIdx.y * blockDim.y + threadIdx.x;
    y = blockIdx.x * blockDim.x + threadIdx.y;
    if (x < n && y < n)
        dst[y*n+x] = data[threadIdx.x][threadIdx.y];
}


/**
 * This function initializes two columns to zero, as necessary for
 * the DST-I
 */
__global__ void init_columns(int n, int m, real *x)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= n) return;

    x[idx*m] = 0.0f;
    x[idx*m+n+1] = 0.0f;
}


/**
 * This function transposes the src vector and stores it
 * in the dst vector.
 * The signal is extended for odd-symmetry. Applying a DFT
 * to this data is necessary to compute the DST-I.
 */
__global__ void copy_flip(int n, int m, real *src, real *dst)
{
    __shared__ real data[BLOCK_SIZE][BLOCK_SIZE+1];

    // column number
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    // row number
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    // read data into shared mem, row-wise (fully coalesced)
    if (x < n && y < n)
        data[threadIdx.y][threadIdx.x] = src[y*n+x];
    
    __syncthreads();

    // copy matrix out-of-place into the buffer
    // and extend signal for odd-symmetry
    if (x  < n && y < n) {
        dst[y*m+1+x] = data[threadIdx.y][threadIdx.x];
        dst[y*m+m-1-x] = -data[threadIdx.y][threadIdx.x];
    }
}

/**
 * Data src is divided by the eigenvalues of the discrete Poisson
 * operator. After that, the data is transposed and stored in
 * dst.
 */
__global__ void scale_copy_flip(int n, int m, real *src, real *dst)
{
    __shared__ real data[BLOCK_SIZE][BLOCK_SIZE];

    // column number
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    // row number
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    // read data to shared mem and scale
    if (x < n && y < n)
        data[threadIdx.y][threadIdx.x] = src[y*n+x]/(4.0 - 2.0*cos(((x+1)*PI)/(n+1)) - 2.0*cos(((y+1)*PI)/(n+1)));

    __syncthreads();

    // copy to global and prepare for FFT
    if (x < n && y < n) {
        dst[y*m+1+x] = data[threadIdx.y][threadIdx.x];
        dst[y*m+m-1-x] = -data[threadIdx.y][threadIdx.x];
    }
}

/**
 * Data src is transposed and scaled by scal as necessary for the DFT.
 * Result is stored in dst.
 */
__global__ void tranpose_scale(int n, int m, real scal, complex *src, real *dst) 
{
    // one column more to avoid bank conflicts
    __shared__ real data[BLOCK_SIZE][BLOCK_SIZE+1];

    // column number
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    // row number
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    // reading imaginary part into shared mem
    if (x < n && y < n)
        data[threadIdx.y][threadIdx.x] = src[y*m+1+x].y*scal;
    
    __syncthreads();

    // transposing matrix
    x = blockIdx.y * blockDim.y + threadIdx.x;
    y = blockIdx.x * blockDim.x + threadIdx.y;
    if (x < n && y < n)
        dst[y*n+x] = data[threadIdx.x][threadIdx.y];
}


/**
 * Extract phase of the realFFT algorithm, which is being used
 * to compute the DFT (of the real signal).
 */
__global__ void extract_scale(int n, int m, real scale, complex *src, real *dst)
{
    __shared__ real data[BLOCK_SIZE][BLOCK_SIZE+1];
    
    // column number
    int x = blockIdx.x * blockDim.x + threadIdx.x + 1;
    // row number
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    
    // global index (for complex type)
    int idx_l = y*(n+1)+x;
    int idx_r = (y+1)*(n+1)-x;

    // read part of data to shared mem
    if (x < n+1 && y < n)
        data[threadIdx.y][threadIdx.x] = src[idx_l].y;

    __syncthreads();

    // extract RFFT to shared mem
    if (x < n+1 && y < n) {
        data[threadIdx.y][threadIdx.x] -= 0.5*((data[threadIdx.y][threadIdx.x]+src[idx_r].y)*(1.0 + sin((2*PI*x)/m)) + (src[idx_l].x-src[idx_r].x)*cos((2*PI*x)/m));
    }
    __syncthreads();

    // transposing matrix
    x = blockIdx.y * blockDim.y + threadIdx.x;
    y = blockIdx.x * blockDim.x + threadIdx.y;
    if (x < n && y < n)
        dst[y*n+x] = data[threadIdx.x][threadIdx.y]*scale;
}

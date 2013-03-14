/* 
 *  main.c - this file is part of CuPoisson
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

/* Driver file as an example for
 * the usage of the poisson solver
 * in your own code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "precision.h"
#include "utils.h"
#include "poisson.h"

int main(int argc, char **argv)
{
    // device memory
    real *psi_d, *z_d;

    size_t fSize = sizeof(real);

    /* grid dimensions */
    unsigned int Nx = 513, Ny = 513;
    // omitting boundaries
    unsigned int nGridPoints = (Nx-2)*(Ny-2);

    cudaMalloc((void **) &psi_d, (nGridPoints+1)*fSize);
    cudaMalloc((void **) &z_d,   (nGridPoints+1)*fSize);

    /* initialization */
    fillArray(psi_d, 0.0, nGridPoints+1);
    fillArray(z_d,   1.0, nGridPoints+1);
    checkCudaError("Initialization of grid");

    // for timing purposes
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    // start timer
    cudaEventRecord(start,0);

    /* Call the poisson solver, right hand side
     * is stored on the device in z_d (make sure the data
     * is copied from CPU to GPU!), result is stored in
     * psi_d (on the GPU/device).
     * Here NX-2 is the width of the grid's interior
     * (without the boundaries).
     */
    cuPoisson((Nx-2), psi_d, z_d);

    // stop timer
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    float computationTime;
    cudaEventElapsedTime(&computationTime, start, stop);

    printf("Computation time was %.5f seconds.\n\n", computationTime/1000.0);

    printf("Writing result to disk...\n");
    // write result to file
    writeBinaryFile(Nx, Ny, psi_d, "data.dat");
    printf("done\n");

    return EXIT_SUCCESS;
}


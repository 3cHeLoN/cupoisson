/* 
 *  poisson.cuh - this file is part of CuPoisson
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


#ifndef POISSON_CUH_

#include "poisson.h"
#include "precision.h"

__global__ void extract_rfft(int n, int m, complex *src, real *dst);
__global__ void init_columns(int n, int m, real *x);
__global__ void copy_flip(int n, int m, real *src, real *dst);
__global__ void scale_copy_flip(int n, int m, real *src, real *dst);
__global__ void tranpose_scale(int n, int m, real scal, complex *src, real *dst);
__global__ void extract_scale(int n, int m, real scale, complex *src, real *dst);

#endif

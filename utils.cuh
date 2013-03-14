/* 
 *  utils.cuh - this file is part of CuPoisson
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


#ifndef UTILS_CUH_
#define UTILS_CUH_

#include "utils.h"

__global__ void cuFillArray(unsigned int n, real *dest, const real value);

#endif

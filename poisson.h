/* 
 *  poisson.h - this file is part of CuPoisson
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

// Include this file in your C source code


#ifndef POISSON_H_
#define POISSON_H_

#include "precision.h"

int cuPoisson(int n, real *u, real *z);
cufftResult cufftExec(cufftHandle plan, complex *idata, complex *odata, int direction);

#endif

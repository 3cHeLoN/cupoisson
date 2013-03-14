#ifndef PRECISION_H_
#define PRECISION_H_

#include <cufft.h>

//#define DOUBLE_PRECISION

#ifdef DOUBLE_PRECISION
    typedef double real;
    typedef cufftDoubleComplex complex;
#else
    typedef float real;
    typedef cufftComplex complex;
#endif

#endif

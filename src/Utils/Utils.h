#ifndef __UTILS_H
#define __UTILS_H

#include "mkl.h"
#include <cmath>

namespace Utils
{
  MKL_INT binsearch(const MKL_INT *array, 
                    MKL_INT len, 
                    MKL_INT value);
  double expectation_value_dense_diag(double *bra,
                                      double *ket,
                                      double *op,
                                      double *tmp,
                                      MKL_INT &basis_size);  
}
#endif

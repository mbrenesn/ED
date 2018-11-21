#ifndef __BASIS_H
#define __BASIS_H

#include "mkl.h"

#include <iostream>
#include <vector>
#include <complex>
#include <math.h>
#include <sys/time.h>

class Basis
{
  public:
    // Methods
    Basis(MKL_INT sites, MKL_INT nup);
    ~Basis();
    Basis(const Basis &rhs);
    Basis &operator=(const Basis &rhs);
    void print_basis();
    // Members
    MKL_INT l, n;
    MKL_INT basis_size;
    MKL_INT *int_basis;
  
  private:
    MKL_INT first_int_();
    MKL_INT basis_size_();
    void construct_int_basis_();
};
#endif

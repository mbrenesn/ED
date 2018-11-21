#ifndef __BASIS_H
#define __BASIS_H

#include <iostream>
#include <vector>
#include <complex>
#include <math.h>
#include <sys/time.h>

class Basis
{
  public:
    // Methods
    Basis(int sites, int nup);
    ~Basis();
    Basis(const Basis &rhs);
    Basis &operator=(const Basis &rhs);
    void print_basis();
    // Members
    int l, n;
    int basis_size;
    int *int_basis;
  
  private:
    int first_int_();
    int basis_size_();
    void construct_int_basis_();
};
#endif

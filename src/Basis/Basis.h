#ifndef __BASIS_H
#define __BASIS_H

#include <iostream>
#include <vector>
#include <complex>
#include <math.h>
#include <sys/time.h>

typedef long long int LLInt;

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
    LLInt basis_size;
    LLInt *int_basis;
  
  private:
    LLInt first_int_();
    LLInt basis_size_();
    void construct_int_basis_();
};
#endif

#include "Basis.h"

Basis::Basis(MKL_INT sites, MKL_INT nup)
: l(sites), n(nup)
{
  basis_size = basis_size_();
  int_basis = new MKL_INT[basis_size];
  construct_int_basis_();
}

Basis::Basis(const Basis &rhs)
{
  std::cout << "Copy constructor (basis) has been called!" << std::endl;

  l = rhs.l;
  n = rhs.n;
  basis_size = rhs.basis_size;
  int_basis = new MKL_INT[basis_size];
  for(MKL_INT i = 0; i < basis_size; ++i)
    int_basis[i] = rhs.int_basis[i];
}

Basis &Basis::operator=(const Basis &rhs)
{
  if(this != &rhs){
    std::cout << "Assignment operator (basis) has been called!" << std::endl;
    delete [] int_basis;
    l = rhs.l;
    n = rhs.n;
    basis_size = rhs.basis_size;
    int_basis = new MKL_INT[basis_size];
    for(MKL_INT i = 0; i < basis_size; ++i)
      int_basis[i] = rhs.int_basis[i];
  }
  
  return *this;
}

Basis::~Basis()
{
  delete [] int_basis;
}

MKL_INT Basis::basis_size_()
{
  double size = 1.0;
  for(MKL_INT i = 1; i <= (l - n); ++i){
    size *= (static_cast<double> (i + n) / static_cast<double> (i));
  }

  return floor(size + 0.5);
}

MKL_INT Basis::first_int_()
{
  MKL_INT first = 0;
  for(MKL_INT i = 0; i < n; ++i){
    first += 1 << i;
  }

  return first;
}

void Basis::construct_int_basis_()
{
  MKL_INT w;
  MKL_INT first = first_int_();

  int_basis[0] = first;

  for(MKL_INT i = 1; i < basis_size; ++i){
    MKL_INT t = (first | (first - 1)) + 1;
    w = t | ((((t & -t) / (first & -first)) >> 1) - 1);
    
    int_basis[i] = w;

    first = w;
  }
}

void Basis::print_basis()
{
  for(MKL_INT i = 0; i <  basis_size; ++i)
    std::cout << int_basis[i] << std::endl;
}

// Intended to be used for calculations on the ETH and Fischer information project
// Computes susceptibility at a given temperature from canonical state,
// the idea is to compare that with the one obtained from f0(E,w) i.e. offdiagonals ETH

#include <iostream>
#include <numeric>
#include <algorithm>
#include <sys/time.h>

#include "mkl.h"
#include "mkl_cblas.h"
#include "mkl_lapacke.h"

#include "../Basis/Basis.h"
#include "../Hamiltonian/XXZ.h"

std::vector<double> linspace(double a, double b, size_t n) {
  double h = (b - a) / static_cast<double>(n - 1);
  std::vector<double> vec(n);
  typename std::vector<double>::iterator x;
  double val;
  for (x = vec.begin(), val = a; x != vec.end(); ++x, val += h)
    *x = val;
  return vec;
}

double seconds()
{
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}

int main(int argc, char **argv)
{
  MKL_INT l = 777;
  MKL_INT n = 777;
  double alpha_val = 0.777;
  double delta_val = 0.777;
  double T = 0.777;
  double eps = 0.777;
  double stag = 0.777;
  bool periodic = true;

  if(argc != 15){
    std::cerr << "Usage: " << argv[0] << 
      " --l [sites] --n [fill] --alpha [alpha] --delta [delta] --eps [eps] --stag [stag] --periodic [bool]" 
        << std::endl;
    exit(1);
  }

  for(int i = 0; i < argc; ++i){
    std::string str = argv[i];
    if(str == "--l") l = atoi(argv[i + 1]);
    else if(str == "--n") n = atoi(argv[i + 1]);
    else if(str == "--alpha") alpha_val = atof(argv[i + 1]);
    else if(str == "--delta") delta_val = atof(argv[i + 1]);
    else if(str == "--eps") eps = atof(argv[i + 1]);
    else if(str == "--stag") stag = atof(argv[i + 1]);
    else if(str == "--periodic") periodic = atoi(argv[i + 1]);
    else continue;
  }

  if(l == 777 || n == 777 || alpha_val == 0.777 || delta_val == 0.777 || eps == 0.777 || stag == 0.777){
    std::cerr << "Error setting parameters" << std::endl;
    std::cerr << "Usage: " << argv[0] << 
      " --l [sites] --n [fill] --alpha [alpha] --delta [delta] --eps [eps] --stag [stag] --periodic [bool]" 
        << std::endl;
    exit(1);
  }

  std::vector<double> alpha(l, alpha_val);
  std::vector<double> delta(l, delta_val);
  std::vector<double> h(l, 0.0);
  // Small perturbation to rid ourselves of symmetries
  h[0] += eps;
  // Staggered field
  for(MKL_INT i = 1; i < l; i += 2) h[i] += -1.0 * stag;

  std::cout << std::fixed;
  std::cout.precision(2);
  std::cout << "# Parameters:" << std::endl;
  std::cout << "# L = " << l << std::endl;
  std::cout << "# N = " << n << std::endl;
  std::cout << "# Periodic boundaries: " << periodic << std::endl;
  std::cout << "# Periodic = 0 means open boundaries, periodic > 0 means periodic boundaries" << std::endl;
  std::cout << "# Alpha = " << "[";
  for(int i = 0; i < (l - 1); ++i){
    std::cout << alpha[i] << ", ";
  }
  std::cout << alpha[l - 1] << "]" << std::endl;
  std::cout << "# Delta = " << "[";
  for(int i = 0; i < (l - 1); ++i){
    std::cout << delta[i] << ", ";
  }
  std::cout << delta[l - 1] << "]" << std::endl;
  std::cout << "# h = " << "[";
  for(int i = 0; i < (l - 1); ++i){
    std::cout << h[i] << ", ";
  }
  std::cout << h[l - 1] << "]" << std::endl;

  Basis *basis = new Basis(l, n);
  //basis->print_basis();

  MKL_INT basis_size = basis->basis_size;

  XXZ heisen( *basis, 
              periodic, 
              true,
              false,
              false );
  heisen.construct_xxz( basis->int_basis, 
                        alpha, 
                        delta, 
                        h);
  //heisen.print_ham_mat();

  delete basis;

  std::cout << std::fixed;
  std::cout.precision(8);
  // Diag
  // TODO Compare this against LAPACKE_dsyevr
  std::vector<double> eigvals(basis_size, 0.0);
  MKL_INT info;

  std::cout << "# Calculating eigen..." << std::endl;
  double tic = seconds();
  info = LAPACKE_dsyevd( LAPACK_ROW_MAJOR,
                         'V', 
                         'U',
                         basis_size, 
                         &heisen.HamMat[0], 
                         basis_size, 
                         &eigvals[0] );
  if(info > 0){
    std::cout << "Eigenproblem" << std::endl;
    exit(1);
  }
  double toc = seconds();
  std::cout << "# Time eigen: " << (toc - tic) << std::endl;


  return 0;
}

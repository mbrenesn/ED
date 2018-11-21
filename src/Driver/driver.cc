#include <iostream>
#include <sys/time.h>

#include "mkl.h"
#include "mkl_lapacke.h"

#include "../Basis/Basis.h"
#include "../Hamiltonian/XXZ.h"

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
  double h_val = 0.777;
  bool periodic = true;

  if(argc != 13){
    std::cerr << "Usage: " << argv[0] << 
      " --l [sites] --n [fill] --alpha [alpha] --delta [delta] --h [h] --periodic [bool]" 
        << std::endl;
    exit(1);
  }

  for(int i = 0; i < argc; ++i){
    std::string str = argv[i];
    if(str == "--l") l = atoi(argv[i + 1]);
    else if(str == "--n") n = atoi(argv[i + 1]);
    else if(str == "--alpha") alpha_val = atof(argv[i + 1]);
    else if(str == "--delta") delta_val = atof(argv[i + 1]);
    else if(str == "--h") h_val = atof(argv[i + 1]);
    else if(str == "--periodic") periodic = atoi(argv[i + 1]);
    else continue;
  }

  if(l == 777 || n == 777 || alpha_val == 0.777 || delta_val == 0.777 || h_val == 0.777){
    std::cerr << "Error setting parameters" << std::endl;
    std::cerr << "Usage: " << argv[0] << 
      " --l [sites] --n [fill] --alpha [alpha] --delta [delta] --h [h] --periodic [bool]" 
        << std::endl;
    exit(1);
  }

  std::vector<double> alpha(l, alpha_val);
  std::vector<double> delta(l, delta_val);
  std::vector<double> h(l, 0.0);
  // Impurity model
  h[l / 2] = h_val;

#if 0
  std::cout << std::fixed;
  std::cout.precision(1);
  std::cout << "# Parameters:" << std::endl;
  std::cout << "# L = " << l << std::endl;
  std::cout << "# N = " << n << std::endl;
  std::cout << "# Periodic boundaries: " << periodic << std::endl;
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
#endif

  Basis *basis = new Basis(l, n);
  //basis->print_basis();

  MKL_INT basis_size = basis->basis_size;

  XXZ heisen( *basis, 
              periodic, 
              true );
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

  //std::cout << "# Calculating eigen..." << std::endl;
  double tic = seconds();
  info = LAPACKE_dsyevd( LAPACK_ROW_MAJOR,
                         'V', 
                         'U',
                         basis_size, 
                         &heisen.HamMat[0], 
                         basis_size, 
                         &eigvals[0] );
  if(info > 0){
    //std::cout << "Eigenproblem" << std::endl;
    exit(1);
  }
  double toc = seconds();
  //std::cout << "# Time eigen: " << (toc - tic) << std::endl;

  // Compute expectation value of SigmaZ[l / 2]
  // Transpose eigvectors
  mkl_dimatcopy( 'R',
                 'T',
                 basis_size,
                 basis_size,
                 1,
                 &heisen.HamMat[0],
                 basis_size,
                 basis_size);

  std::vector<double> tmp(basis_size, 0.0);
  std::vector<double> mags(basis_size, 0.0);
  for(MKL_INT i = 0; i < basis_size; ++i){
    vdMul(basis_size,
          &heisen.HamMat[(i * basis_size)],
          &heisen.SigmaZ[l / 2][0],
          &tmp[0]);
    double val = cblas_ddot( basis_size,
                             &heisen.HamMat[(i * basis_size)],
                             1,
                             &tmp[0],
                             1);
    std::cout << (eigvals[i] / l) << " " << val << std::endl;
  }

  return 0;
}

// Impurity calculations ETH

#include <iostream>
#include <sys/time.h>

#include "mkl.h"
#include "mkl_lapacke.h"
#include "mkl_spblas.h"

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
  double eps = 0.777;
  double stag = 0.777;
  bool periodic = true;

  if(argc != 17){
    std::cerr << "Usage: " << argv[0] << 
      " --l [sites] --n [fill] --alpha [alpha] --delta [delta] --h [h] --eps [eps] --stag [stag] --periodic [bool]" 
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
    else if(str == "--eps") eps = atof(argv[i + 1]);
    else if(str == "--stag") stag = atof(argv[i + 1]);
    else if(str == "--periodic") periodic = atoi(argv[i + 1]);
    else continue;
  }

  if(l == 777 || n == 777 || alpha_val == 0.777 || delta_val == 0.777 || h_val == 0.777 || eps == 0.777 || stag == 0.777){
    std::cerr << "Error setting parameters" << std::endl;
    std::cerr << "Usage: " << argv[0] << 
      " --l [sites] --n [fill] --alpha [alpha] --delta [delta] --h [h] --eps [eps] --stag [stag] --periodic [bool]" 
        << std::endl;
    exit(1);
  }

  std::vector<double> alpha(l, alpha_val);
  std::vector<double> delta(l, delta_val);
  std::vector<double> h(l, 0.0);
  // Small perturbation to rid ourselves of symmetries
  h[0] += eps;
  // Impurity model
  h[(l / 2) - 1] += h_val;
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

  // Compute expectation values
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
  // Diagonal matrix for \sigma^z(N/2) * \sigma^z(N/2 + 1)
  std::vector<double> Sn2Sn21(basis_size, 0.0);
  vdMul(basis_size,
        &heisen.SigmaZ[(l / 2) - 1][0],
        &heisen.SigmaZ[(l / 2)][0],
        &Sn2Sn21[0]);
  // Diagonal matrix for \sigma^z(N/4) * \sigma^z(N/4 + 1)
  std::vector<double> Sn4Sn41(basis_size, 0.0);
  vdMul(basis_size,
        &heisen.SigmaZ[(l / 4) - 1][0],
        &heisen.SigmaZ[(l / 4)][0],
        &Sn4Sn41[0]);
  // Diagonal matrix for \sum_i \sigma^z(i) * \sigma^z(i + 1)
  std::vector<double> nn_corr(basis_size, 0.0);
  for(MKL_INT i = 0; i < basis_size; ++i){
    double l_corr = 0.0;
    for(MKL_INT sp = 0; sp < (l - 1); ++sp){
      l_corr += heisen.SigmaZ[sp][i] * heisen.SigmaZ[sp + 1][i];
    }
    nn_corr[i] = l_corr / (l - 1);
  }

  // Off diagonals
  double window = 0.05;
  
  std::cout << std::scientific;
  tic = seconds();
  for(MKL_INT i = 0; i < basis_size; ++i){
    for(MKL_INT j = 0; j < (i + 1); ++j){
      if( i == j ) continue;
      if( ( (std::abs(eigvals[i] + eigvals[j])) / l ) <= window ){
        // Sigma^z_(N/4) * Sigma^z_(N/4 + 1)
        double val1 = Utils::expectation_value_dense_diag( &heisen.HamMat[(i * basis_size)],
                                                           &heisen.HamMat[(j * basis_size)],
                                                           &Sn4Sn41[0],
                                                           &tmp[0],
                                                           basis_size);
        // Sigma^z_(N/2) * Sigma^z_(N/2 + 1)
        double val2 = Utils::expectation_value_dense_diag( &heisen.HamMat[(i * basis_size)],
                                                           &heisen.HamMat[(j * basis_size)],
                                                           &Sn2Sn21[0],
                                                           &tmp[0],
                                                           basis_size);
        // \sum_i^{L-1} Sigma^z_(i) * Sigma^z_(i + 1)
        double val3 = Utils::expectation_value_dense_diag( &heisen.HamMat[(i * basis_size)],
                                                           &heisen.HamMat[(j * basis_size)],
                                                           &nn_corr[0],
                                                           &tmp[0],
                                                           basis_size);
        std::cout << eigvals[i] - eigvals[j] << " " << std::abs(val1) << " " 
          << std::abs(val2) << " " << std::abs(val3) << std::endl;
      }
    }
  }
  toc = seconds();
  std::cout << "# Time mult: " << (toc - tic) << std::endl;

  return 0;
}

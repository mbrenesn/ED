// Stark localisation

#include <iostream>
#include <sys/time.h>

#include "mkl.h"
#include "mkl_lapacke.h"

#include "../Basis/Basis.h"
#include "../Hamiltonian/StarkM.h"

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
  double J = 0.777;
  double V = 0.777;
  double gamma = 0.777;
  double alpha = 0.777;

  if(argc != 13){
    std::cerr << "Usage: " << argv[0] << 
      " --l [sites] --n [fill] --J [J] --V [V] --gamma [gamma] --alpha [alpha]" 
        << std::endl;
    exit(1);
  }

  for(int i = 0; i < argc; ++i){
    std::string str = argv[i];
    if(str == "--l") l = atoi(argv[i + 1]);
    else if(str == "--n") n = atoi(argv[i + 1]);
    else if(str == "--J") J = atof(argv[i + 1]);
    else if(str == "--V") V = atof(argv[i + 1]);
    else if(str == "--gamma") gamma = atof(argv[i + 1]);
    else if(str == "--alpha") alpha = atof(argv[i + 1]);
    else continue;
  }

  if(l == 777 || n == 777 || J == 0.777 || V == 0.777 || gamma == 0.777 || alpha == 0.777){
    std::cerr << "Error setting parameters" << std::endl;
    std::cerr << "Usage: " << argv[0] << 
      " --l [sites] --n [fill] --J [J] --V [V] --gamma [gamma] --alpha [alpha]" 
        << std::endl;
    exit(1);
  }

  std::cout << std::fixed;
  std::cout.precision(4);
#if 0  
  std::cout << "# Parameters:" << std::endl;
  std::cout << "# L = " << l << std::endl;
  std::cout << "# N = " << n << std::endl;
  std::cout << "# J = " << J << std::endl;
  std::cout << "# V = " << V << std::endl;
  std::cout << "# gamma = " << gamma << std::endl;
  std::cout << "# alpha = " << alpha << std::endl;
#endif 

  Basis *basis = new Basis(l, n);
  //basis->print_basis();

  MKL_INT basis_size = basis->basis_size;

  StarkM stark( *basis, 
                false );
  stark.construct_starkm( basis->int_basis, 
                          J,
                          V,
                          gamma,
                          alpha);
  stark.print_ham_mat();

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
                         &stark.HamMat[0], 
                         basis_size, 
                         &eigvals[0] );
  if(info > 0){
    std::cout << "Eigenproblem" << std::endl;
    exit(1);
  }
  double toc = seconds();

  for(MKL_INT i = 0; i < basis_size; ++i){
    double ent = 0.0;
    for(MKL_INT j = 0; j < basis_size; ++j){
      double pi = std::pow(stark.HamMat[(j * basis_size) + i], 2.0);
      ent += -1.0 * pi * std::log(pi);
    }
    std::cout << alpha << " " 
      << (eigvals[i] - eigvals[basis_size - 1]) / (eigvals[0] - eigvals[basis_size - 1]) 
        << " " << ent << std::endl;
  }
  
  //std::cout << "# Time eigen: " << (toc - tic) << std::endl;

  return 0;
}

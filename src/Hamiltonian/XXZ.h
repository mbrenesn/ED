/* 
 * XXZ Hamiltonian
 * Currently defined only for spin 1/2 systems, periodic or open boundaries
 * SigmaZ mats implemented
 */

#ifndef __XXZ_H
#define __XXZ_H

#include "../Basis/Basis.h"
#include "../Utils/Utils.h"

class XXZ
{
  public:
  // Methods
  XXZ(Basis &basis,
      bool periodic = false,
      bool sigma_z_mats = false,
      bool local_kinetic = false);
  ~XXZ();
  void construct_xxz(MKL_INT *int_basis,
                     std::vector<double> &alpha,
                     std::vector<double> &delta,
                     std::vector<double> &h);
  void print_ham_mat();
  // Members
  std::vector<double> HamMat;
  std::vector< std::vector<double> > SigmaZ;
  std::vector<double> LocalK_vals;
  std::vector<MKL_INT> LocalK_cols;
  std::vector<MKL_INT> LocalK_rowptr;

  private:
    MKL_INT l_; 
    MKL_INT n_;
    bool periodic_;
    bool sigma_z_mats_;
    bool local_kinetic_;
    MKL_INT basis_size_;
};

#endif

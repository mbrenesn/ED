/* 
 * Stark localisation model
 * Currently defined only for spinless fermions and open boundaries
 * Fermion occupation mats implemented
 * Hamiltonian from Eq 1 in PRL 122, 040606 (2019)
 */

#ifndef __STARKM_H
#define __STARKM_H

#include "../Basis/Basis.h"
#include "../Utils/Utils.h"

class StarkM
{
  public:
  // Methods
  StarkM(Basis &basis,
      bool occupation_mats = false);
  ~StarkM();
  void construct_starkm(MKL_INT *int_basis,
                        double &J,
                        double &V,
                        double &gamma,
                        double &alpha);
  void print_ham_mat();
  // Members
  std::vector<double> HamMat;
  std::vector< std::vector<double> > NumF;

  private:
    MKL_INT l_; 
    MKL_INT n_;
    bool occupation_mats_;
    MKL_INT basis_size_;
};
#endif

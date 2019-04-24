#include "StarkM.h"

StarkM::StarkM(Basis &basis,
               bool occupation_mats)
: l_(basis.l), n_(basis.n), basis_size_(basis.basis_size), occupation_mats_(occupation_mats)
{
  HamMat.resize(basis_size_ * basis_size_);
  if(occupation_mats_){
    NumF.resize(l_);
    for(MKL_INT i = 0; i < l_; ++i)
      NumF[i].resize(basis_size_);
  }
}

StarkM::~StarkM()
{}

void StarkM::construct_starkm(MKL_INT *int_basis,
                              double &J,
                              double &V,
                              double &gamma,
                              double &alpha)
{
  for(MKL_INT state = 0; state < basis_size_; ++state){

    MKL_INT bs = int_basis[state];
    double vi = 0.0;
    double local_term = 0.0;

    for(MKL_INT site = 0; site < l_; ++site){
      MKL_INT bitset = bs;
      if(bitset & (1LL << site)){
        if(occupation_mats_)
          NumF[site][state] = 1.0;
        local_term += ((-1.0 * gamma * site) + (alpha * site * site / (l_ * l_))) * 0.5;
      }
      else{
        if(occupation_mats_)
          NumF[site][state] = 0.0;
        local_term -= ((-1.0 * gamma * site) + (alpha * site * site / (l_ * l_))) * 0.5;
      }

      if(site == l_ - 1) continue;

      if(bitset & (1LL << site)){
        if(bitset & (1LL << ( (site + 1) % l_ ) ) ){
          vi += V * 0.25;
          continue;
        }
        else{
          vi -= V * 0.25;
          bitset ^= 1LL << site;
          bitset ^= 1LL << (site + 1) % l_;
          MKL_INT match_ind1 = Utils::binsearch(int_basis, basis_size_, bitset);
          HamMat[ (state * basis_size_) + match_ind1 ] = 0.5 * J;
          continue;
        }     
      } // End spin up case
      else{
        if(bitset & (1LL << ( (site + 1) % l_ ) ) ){
          vi -= V * 0.25;
          bitset ^= 1LL << site;
          bitset ^= 1LL << (site + 1) % l_;
          MKL_INT match_ind2 = Utils::binsearch(int_basis, basis_size_, bitset);
          HamMat[ (state * basis_size_) + match_ind2 ] = 0.5 * J;
          continue;
        }
        else{
          vi += V * 0.25;
          continue;
        }
      } // End spin down case
    } // End site loop
  HamMat[ (state * basis_size_) + state ] = vi + local_term;
  } // End state loop
}

void StarkM::print_ham_mat()
{
  for(MKL_INT i = 0; i < basis_size_; ++i){
    for(MKL_INT j = 0; j < basis_size_; ++j){
      std::cout << HamMat[(i * basis_size_) + j] << " ";
    }
    std::cout << std::endl; 
  }
}

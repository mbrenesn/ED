#include "XXZ.h"

XXZ::XXZ(Basis &basis,
         bool periodic,
         bool sigma_z_mats)
: l_(basis.l), n_(basis.n), basis_size_(basis.basis_size), 
periodic_(periodic), sigma_z_mats_(sigma_z_mats)
{
  HamMat.resize(basis_size_ * basis_size_);
  if(sigma_z_mats_){
    SigmaZ.resize(l_);
    for(int i = 0; i < l_; ++i)
      SigmaZ[i].resize(basis_size_);
  }
}

XXZ::~XXZ()
{}

void XXZ::construct_xxz(LLInt *int_basis,
                        std::vector<double> &alpha,
                        std::vector<double> &delta,
                        std::vector<double> &h)
{
  for(int state = 0; state < basis_size_; ++state){

    LLInt bs = int_basis[state];
    double vi = 0.0;
    double mag_term = 0.0;

    for(int site = 0; site < l_; ++site){
      int bitset = bs;
      if(bitset & (1 << site)){
        if(sigma_z_mats_)
          SigmaZ[site][state] = 1.0;
        mag_term += h[site];
      }
      else{
        if(sigma_z_mats_)
          SigmaZ[site][state] = -1.0;
        mag_term -= h[site];
      }

      if(!periodic_){
        if(site == l_ - 1)
          continue;
      }

      if(bitset & (1 << site)){
        if(bitset & (1 << ( (site + 1) % l_ ) ) ){
          vi += delta[site];
          continue;
        }
        else{
          vi -= delta[site];
          bitset ^= 1 << site;
          bitset ^= 1 << (site + 1) % l_;
          LLInt match_ind1 = Utils::binsearch(int_basis, basis_size_, bitset);
          HamMat[ (state * basis_size_) + match_ind1 ] = 2.0 * alpha[site];
          continue;
        }     
      } // End spin up case
      else{
        if(bitset & (1 << ( (site + 1) % l_ ) ) ){
          vi -= delta[site];
          bitset ^= 1 << site;
          bitset ^= 1 << (site + 1) % l_;
          LLInt match_ind2 = Utils::binsearch(int_basis, basis_size_, bitset);
          HamMat[ (state * basis_size_) + match_ind2 ] = 2.0 * alpha[site];
          continue;
        }
        else{
          vi += delta[site];
          continue;
        }
      } // End spin down case
    } // End site loop
  HamMat[ (state * basis_size_) + state ] = vi + mag_term;
  } // End state loop
}

void XXZ::print_ham_mat()
{
  for(LLInt i = 0; i < basis_size_; ++i){
    for(LLInt j = 0; j < basis_size_; ++j){
      std::cout << HamMat[(i * basis_size_) + j] << " ";
    }
    std::cout << std::endl; 
  }
}

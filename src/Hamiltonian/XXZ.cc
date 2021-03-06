#include "XXZ.h"

XXZ::XXZ(Basis &basis,
         bool periodic,
         bool sigma_z_mats,
         bool local_kinetic,
         bool current)
: l_(basis.l), n_(basis.n), basis_size_(basis.basis_size), 
periodic_(periodic), sigma_z_mats_(sigma_z_mats), local_kinetic_(local_kinetic), current_(current)
{
  HamMat.resize(basis_size_ * basis_size_);
  if(sigma_z_mats_){
    SigmaZ.resize(l_);
    for(MKL_INT i = 0; i < l_; ++i)
      SigmaZ[i].resize(basis_size_);
  }
  if(local_kinetic_)
    LocalK_rowptr.push_back(0);
  if(current_)
    Curr_rowptr.push_back(0);
}

XXZ::~XXZ()
{}

void XXZ::construct_xxz(MKL_INT *int_basis,
                        std::vector<double> &alpha,
                        std::vector<double> &delta,
                        std::vector<double> &h)
{
  for(MKL_INT state = 0; state < basis_size_; ++state){

    MKL_INT bs = int_basis[state];
    double vi = 0.0;
    double mag_term = 0.0;

    for(MKL_INT site = 0; site < l_; ++site){
      MKL_INT bitset = bs;
      if(bitset & (1LL << site)){
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

      if(bitset & (1LL << site)){
        if(bitset & (1LL << ( (site + 1) % l_ ) ) ){
          vi += delta[site];
          continue;
        }
        else{
          vi -= delta[site];
          bitset ^= 1LL << site;
          bitset ^= 1LL << (site + 1) % l_;
          MKL_INT match_ind1 = Utils::binsearch(int_basis, basis_size_, bitset);
          HamMat[ (state * basis_size_) + match_ind1 ] = 2.0 * alpha[site];
          if( local_kinetic_ && ( site == ((l_ / 2) - 1) ) ){ 
            LocalK_vals.push_back( 2.0 * alpha[site] );
            LocalK_cols.push_back( match_ind1 );
          }
          if( current_ ){
            Curr_vals.push_back( -2.0 * alpha[site] );
            Curr_cols.push_back( match_ind1 );
          }
          continue;
        }     
      } // End spin up case
      else{
        if(bitset & (1LL << ( (site + 1) % l_ ) ) ){
          vi -= delta[site];
          bitset ^= 1LL << site;
          bitset ^= 1LL << (site + 1) % l_;
          MKL_INT match_ind2 = Utils::binsearch(int_basis, basis_size_, bitset);
          HamMat[ (state * basis_size_) + match_ind2 ] = 2.0 * alpha[site];
          if( local_kinetic_ && ( site == ((l_ / 2) - 1) ) ){ 
            LocalK_vals.push_back( 2.0 * alpha[site] );
            LocalK_cols.push_back( match_ind2 );
          }
          if( current_ ){
            Curr_vals.push_back( 2.0 * alpha[site] );
            Curr_cols.push_back( match_ind2 );
          }
          continue;
        }
        else{
          vi += delta[site];
          continue;
        }
      } // End spin down case
    } // End site loop
    HamMat[ (state * basis_size_) + state ] = vi + mag_term;
    if(local_kinetic_) LocalK_rowptr.push_back( LocalK_vals.size() );
    if(current_) Curr_rowptr.push_back( Curr_vals.size() );
  } // End state loop
}

void XXZ::print_ham_mat()
{
  for(MKL_INT i = 0; i < basis_size_; ++i){
    for(MKL_INT j = 0; j < basis_size_; ++j){
      std::cout << HamMat[(i * basis_size_) + j] << " ";
    }
    std::cout << std::endl; 
  }
}

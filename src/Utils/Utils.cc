#include "Utils.h"

namespace Utils
{
  /*******************************************************************************/
  // Binary search: Divide and conquer. For the construction of the Hamiltonian
  // matrix instead of looking through all the elements of the int basis a
  // binary search will perform better for large systems
  /*******************************************************************************/
  MKL_INT binsearch(const MKL_INT *array, MKL_INT len, MKL_INT value)
  {
    if(len == 0) return -1;
    MKL_INT mid = len / 2;

    if(array[mid] == value) 
      return mid;
    else if(array[mid] < value){
      MKL_INT result = binsearch(array + mid + 1, len - (mid + 1), value);
      if(result == -1) 
        return -1;
      else
        return result + mid + 1;
    }
    else
      return binsearch(array, mid, value);
  }
}

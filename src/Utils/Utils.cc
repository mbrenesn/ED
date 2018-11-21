#include "Utils.h"

namespace Utils
{
  /*******************************************************************************/
  // Binary search: Divide and conquer. For the construction of the Hamiltonian
  // matrix instead of looking through all the elements of the int basis a
  // binary search will perform better for large systems
  /*******************************************************************************/
  LLInt binsearch(const LLInt *array, LLInt len, LLInt value)
  {
    if(len == 0) return -1;
    LLInt mid = len / 2;

    if(array[mid] == value) 
      return mid;
    else if(array[mid] < value){
      LLInt result = binsearch(array + mid + 1, len - (mid + 1), value);
      if(result == -1) 
        return -1;
      else
        return result + mid + 1;
    }
    else
      return binsearch(array, mid, value);
  }
}

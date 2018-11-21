#ifndef __UTILS_H
#define __UTILS_H

#include <cmath>

typedef long long int LLInt;

namespace Utils
{
  LLInt binsearch(const LLInt *array, 
                  LLInt len, 
                  LLInt value);
}
#endif

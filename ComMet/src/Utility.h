// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#ifndef __INC_UTILITY_H__
#define __INC_UTILITY_H__

#include <utility>
#include <string>

#include "LogSpace.hpp"

#define PI 3.141592653589793238462643383279502884197169399375105820974944

#ifndef uint
#define uint unsigned int
#endif

void
progress(const char* message);

uint 
round(float val);

double
log_binom(uint m, uint u, double s);

double
log_asin(double ds);

double
log_beta(double ds, double alpha);

#endif

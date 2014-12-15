#ifndef __INC_UTILITY_HH__
#define __INC_UTILITY_HH__

#include <string>

#include "LogSpace.hh"

#define PI    3.141592653589793238462643383279502884197169399375105820974944
#define EULER 0.577215664901532860606512090082402431042159335939923598805767
#define RAND_SEED 1985

#ifndef uint
#define uint unsigned int 
#endif

std::string
whoami();

void
progress(const char* message);

void
progress(const std::string& message);

// instead of clumsy <iomanip> ...
std::string 
formatval(const char* format, double val);

uint 
round(float x);

double 
fast_gamma(double x);

double 
fast_digamma(double x);

double 
fast_beta(double x, double y);

double 
fast_comb(double x, double y);

double 
log_binom(uint m, uint u, double s);

#endif

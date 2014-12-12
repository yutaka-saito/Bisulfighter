#include <iostream>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "Utility.hh"

using namespace std;
using namespace boost;
namespace bm = boost::math;

std::string
whoami()
{
  const string iam = "ComMet v1.1 - identification of differentially methylated regions (DMRs)";
  return iam;
}

void
progress(const char* message) 
{
  cout << message << endl << flush;
}

void
progress(const string& message) 
{
  cout << message << endl << flush;
}

string 
formatval(const char* format, double val)
{
  char buf[64];
  sprintf(buf, format, val);
  return string(buf);
}

uint 
round(float x) 
{
  assert(x >= 0.0);
  return (uint) (x + 0.5);
}

double 
fast_gamma(double x) 
{
  assert(x > 0.0);
  return bm::lgamma(x);
}

double 
fast_digamma(double x) 
{
  assert(x > 0.0);
  return bm::digamma(x);
}

double 
fast_beta(double x, double y)
{
  assert(x > 0.0 && y > 0.0);
  return log(bm::beta(x, y));
  //return bm::lgamma(x) + bm::lgamma(x) - bm::lgamma(x+y);
}

double 
fast_comb(double x, double y)
{
  assert(x >= 0.0 && y >= 0.0);
  return log(bm::binomial_coefficient<double>((uint) x, (uint) y));
  //return bm::lgamma(x+1.0) - bm::lgamma(y+1.0) - bm::lgamma(x-y+1.0);
}

double 
log_binom(uint m, uint u, double s) {
  uint n = m + u;
  double ls = log(s);
  double lf = log(1.0 - s);

  if (m==0) { 
    return n * lf;
  }
  else if (u==0) {
    return n * ls;
  }
  else {
    double lprb = m * ls + u * lf;
    for (uint i=0; i!=m; ++i) 
      lprb += log((double) (u + i + 1) / (i + 1));
    return lprb;
  }
}

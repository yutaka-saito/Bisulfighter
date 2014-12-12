// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#include <iostream>
#include <cassert>
#include <cmath>

#include "Utility.h"

using namespace std;

void
progress(const char* message) {
  cout << message << endl << flush;
}

void
progress(const string& message) {
  cout << message << endl << flush;
}

string 
whoami() {
  string iam = "ComMet v1.0 - identification of differentially methylated regions (DMRs)";
  return iam;
}

uint 
round(float val) {
  return (uint) (val + 0.5);
}

double 
log_binom(uint m, uint u, double s) {
  uint n = m + u;
  double ls = Log(s);
  double lf = Log(1.0 - s);

  if (m==0) { 
    return n * lf;
  }
  else if (u==0) {
    return n * ls;
  }
  else {
    double lprb = m * ls + u * lf;
    for (uint i=0; i!=m; ++i) 
      lprb += Log((double) (u + i + 1) / (i + 1));
    return lprb;
  }
}

double
log_asin(double ds) {
  ds = 0.5 + (ds + 1.0) / 4;
  return Log(2.0 / (PI * sqrt(ds * (1.0 - ds) ) ) );
}

double
log_beta(double ds, double alpha) {
  ds = (ds + 1.0) / 2;
  return Log(alpha * pow(ds, alpha - 1.0)) ;
}


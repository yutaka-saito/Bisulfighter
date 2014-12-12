#ifndef __INC_GLOBAL_STATISTICS_HH__
#define __INC_GLOBAL_STATISTICS_HH__

#include <vector>
#include <cassert>

#include "Utility.hh"
#include "Data.hh"

class GlobalStatistics 
{
public:
  typedef double ValueType;

public:
  GlobalStatistics() {}
  ~GlobalStatistics() {}

public:
  bool reset(const std::vector<MethylList>& data, 
	     uint nmix, bool nobeta, uint nitr, uint nsmp, bool verbose);

  ValueType pup(const MethylList& met, uint i) const;
  ValueType pdown(const MethylList& met, uint i) const;
  ValueType pnochange(const MethylList& met, uint i) const;
  ValueType pup(float m1, float u1, float m2, float u2) const;
  ValueType pdown(float m1, float u1, float m2, float u2) const;
  ValueType pnochange(float m1, float u1, float m2, float u2) const;

  void print_param(uint c);
  void test(bool nobeta);

private:
  void reset_beta_mixture(const std::vector<MethylList>& data, uint c, 
			  uint nitr, uint nsmp, bool verbose);

  ValueType pup_beta_mixture(const MethylList& met, uint i) const;
  ValueType pdown_beta_mixture(const MethylList& met, uint i) const;
  ValueType pnochange_beta_mixture(const MethylList& met, uint i) const;
  ValueType pup_pseudo_count(const MethylList& met, uint i) const;
  ValueType pdown_pseudo_count(const MethylList& met, uint i) const;
  ValueType pnochange_pseudo_count(const MethylList& met, uint i) const;
  ValueType pup_beta_mixture(float m1, float u1, float m2, float u2) const;
  ValueType pdown_beta_mixture(float m1, float u1, float m2, float u2) const;
  ValueType pnochange_beta_mixture(float m1, float u1, float m2, float u2) const;
  ValueType pup_pseudo_count(float m1, float u1, float m2, float u2) const;
  ValueType pdown_pseudo_count(float m1, float u1, float m2, float u2) const;
  ValueType pnochange_pseudo_count(float m1, float u1, float m2, float u2) const;

  ValueType gamma(float av) const;
  ValueType digamma(float av) const;
  ValueType beta(float av, float bv) const;
  ValueType comb(float av, float bv) const; 
  ValueType gamma(uint a) const;
  ValueType digamma(uint a) const;
  ValueType beta(uint a, uint b) const;
  ValueType comb(uint a, uint b) const; 

private:
  uint nmix_;
  bool nobeta_;
  std::vector<ValueType> wh1_;
  std::vector<float> ah1_;
  std::vector<float> bh1_;
  std::vector<ValueType> wl1_;
  std::vector<float> al1_;
  std::vector<float> bl1_;
  std::vector<ValueType> wh2_;
  std::vector<float> ah2_;
  std::vector<float> bh2_;
  std::vector<ValueType> wl2_;
  std::vector<float> al2_;
  std::vector<float> bl2_;
  uint gamma_size_;
  std::vector<ValueType> gamma_;
  uint digamma_size_;
  std::vector<ValueType> digamma_;
};

#endif

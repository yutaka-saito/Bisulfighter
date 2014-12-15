// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#ifndef __INC_NAIVE_MODEL_H__
#define __INC_NAIVE_MODEL_H__

#include <fstream>
#include <vector>

#include "Utility.h"
#include "ProbabilityModel.h"
#include "Data.h"

class NaiveModel : public ProbabilityModel 
{
private:
  enum STATE {
    CPG_UP        =  0,
    CPG_DOWN      =  1,
    CPG_NOCHANGE  =  2,
    GAP_UP        =  3,
    GAP_DOWN      =  4,
    GAP_NOCHANGE  =  5,
    N_STATE       =  6,
    N_CPG_STATE   =  3,
    N_GAP_STATE   =  3,
  };

public:
  NaiveModel() : ProbabilityModel(N_STATE, N_CPG_STATE, N_GAP_STATE) {}
  ~NaiveModel() {}

public:
  void reset_param(const MethylList& met, ValueType alpha);
  void dbase(std::ofstream& ofs, const MethylList& met);
  void dregion(std::ofstream& ofs, const MethylList& met, const std::vector<uint>& path);
  void dregion_viterbi(std::ofstream& ofs, const MethylList& met);
  void dregion_pdecode(std::ofstream& ofs, const MethylList& met);
  void dregion_ratio(std::ofstream& ofs, const MethylList& met, ValueType thsh);
};

#endif	

#ifndef __INC_NAIVE_MODEL_HH__
#define __INC_NAIVE_MODEL_HH__

#include <fstream>
#include <vector>

#include "Utility.hh"
#include "ProbabilityModel.hh"
#include "GlobalStatistics.hh"
#include "Data.hh"

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
  void reset_param(const MethylList& met, const GlobalStatistics& gstat);
  void dbase(std::ofstream& ofs, const MethylList& met);
  void dregion(std::ofstream& ofs, const MethylList& met, const std::vector<uint>& path);
  void dregion_viterbi(std::ofstream& ofs, const MethylList& met);
  void dregion_pdecode(std::ofstream& ofs, const MethylList& met);
  void dregion_ratio(std::ofstream& ofs, const MethylList& met, ValueType thsh);
};

#endif	

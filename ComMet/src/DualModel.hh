#ifndef __INC_DUAL_MODEL_HH__
#define __INC_DUAL_MODEL_HH__

#include <fstream>
#include <vector>

#include "Utility.hh"
#include "ProbabilityModel.hh"
#include "GlobalStatistics.hh"
#include "Data.hh"

class DualModel : public ProbabilityModel 
{
private:
  enum STATE {
    CPG_UP         =  0,
    CPG_DOWN       =  1,
    CPG_NOCHANGE   =  2,
    CPG_UP2        =  3,
    CPG_DOWN2      =  4,
    CPG_NOCHANGE2  =  5,
    GAP_UP         =  6,
    GAP_DOWN       =  7,
    GAP_NOCHANGE   =  8,
    GAP_UP2        =  9,
    GAP_DOWN2      =  10,
    GAP_NOCHANGE2  =  11,
    N_STATE        =  12,
    N_CPG_STATE    =  6,
    N_GAP_STATE    =  6,
  };

public:
  DualModel() : ProbabilityModel(N_STATE, N_CPG_STATE, N_GAP_STATE) {}
  ~DualModel() {}

public:
  void reset_param(const MethylList& met, const GlobalStatistics& gstat);
  void dbase(std::ofstream& ofs, const MethylList& met);
  void dregion(std::ofstream& ofs, const MethylList& met, const std::vector<uint>& path);
  void dregion_viterbi(std::ofstream& ofs, const MethylList& met);
  void dregion_pdecode(std::ofstream& ofs, const MethylList& met);
  void dregion_ratio(std::ofstream& ofs, const MethylList& met, ValueType thsh);
};

#endif	

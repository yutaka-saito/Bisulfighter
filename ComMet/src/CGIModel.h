// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#ifndef __INC_CGI_MODEL_H__
#define __INC_CGI_MODEL_H__

#include <fstream>
#include <vector>

#include "Utility.h"
#include "ProbabilityModel.h"
#include "Data.h"

class CGIModel : public ProbabilityModel 
{
private:
  enum STATE {
    CPG_UP            =  0,
    CPG_DOWN          =  1,
    CPG_NOCHANGE      =  2,
    CGI_CPG_UP        =  3,
    CGI_CPG_DOWN      =  4,
    CGI_CPG_NOCHANGE  =  5,
    GAP_UP            =  6,
    GAP_DOWN          =  7,
    GAP_NOCHANGE      =  8,
    CGI_GAP_UP        =  9,
    CGI_GAP_DOWN      =  10,
    CGI_GAP_NOCHANGE  =  11,
    N_STATE           =  12,
    N_CPG_STATE       =  6,
    N_GAP_STATE       =  6,
  };

public:
  CGIModel() : ProbabilityModel(N_STATE, N_CPG_STATE, N_GAP_STATE) {}
  ~CGIModel() {}

public:
  void reset_param(const MethylList& met, ValueType alpha);
  void dbase(std::ofstream& ofs, const MethylList& met);
  void dregion(std::ofstream& ofs, const MethylList& met, const std::vector<uint>& path);
  void dregion_viterbi(std::ofstream& ofs, const MethylList& met);
  void dregion_pdecode(std::ofstream& ofs, const MethylList& met);
  void dregion_ratio(std::ofstream& ofs, const MethylList& met, ValueType thsh);
};

#endif	

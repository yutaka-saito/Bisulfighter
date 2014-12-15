// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#ifndef __INC_PROBABILITY_MODEL_H__
#define __INC_PROBABILITY_MODEL_H__

#include <fstream>
#include <vector>

#include "Utility.h"
#include "Data.h"

class ProbabilityModel {
protected:
  typedef double ValueType; // float causes numerical errors

protected:
  ProbabilityModel(uint n, uint ncpg, uint ngap) 
    : NState(n), NCPGState(ncpg), NGAPState(ngap) {}
  virtual ~ProbabilityModel() {}

public:
  bool reset(const MethylList& met, ValueType alpha);
  virtual void reset_param(const MethylList& met, ValueType alpha) = 0;
  void print_param(bool logspc);
  void print_table(const std::vector<std::vector<ValueType> >& tbl, bool logspc);

  void viterbi();
  void pdecode(ValueType gcpg, ValueType ggap);
  void trace(std::vector<uint>& path, 
	     const std::vector<std::vector<ValueType> >& tbl,
	     const std::vector<std::vector<uint> >& trc_tbl);
  void trace_viterbi();
  void trace_pdecode();
  void posterior();
  void fwdbwd();
  void em(uint nitr, bool verbose);

protected:
  const uint NState;
  const uint NCPGState;
  const uint NGAPState;
  std::vector<ValueType> InitProb;
  std::vector<std::vector<ValueType> > TransProb;
  std::vector<std::vector<ValueType> > EmitProb;

protected:
  uint dpsize_;
  std::vector<std::vector<ValueType> > vtb_; // DP table for Viterbi algorithm
  std::vector<std::vector<uint> > trc_vtb_; // traceback table for Viterbi algorithm
  std::vector<uint> path_vtb_; // optimal path in trc_vtb_
  std::vector<std::vector<ValueType> > pdc_; // DP table for posterior deconding
  std::vector<std::vector<uint> > trc_pdc_; // traceback table for posterior decoding
  std::vector<uint> path_pdc_; // optimal path in trc_pdc_
  std::vector<std::vector<ValueType> > ppr_; // posterior probability values
  std::vector<std::vector<ValueType> > fwd_; // DP table for forward algorithm
  std::vector<std::vector<ValueType> > bwd_; // DP table for backward algorithm
  ValueType total_;

  bool flg_vtb_;
  bool flg_trc_vtb_;
  bool flg_pdc_;
  bool flg_trc_pdc_;
  bool flg_ppr_;
  bool flg_fwd_;
  bool flg_bwd_;
};

#endif	

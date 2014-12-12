// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#include <iostream>
#include <cassert>

#include "ProbabilityModel.h"

using namespace std;

bool ProbabilityModel::
reset(const MethylList& met, ValueType alpha)
{
  // dpsize_
  dpsize_ = met.pos_.back() - met.pos_[0] + 1;
  cout << "DP size: " << dpsize_ << endl;
  if (dpsize_ < 2 || met.pos_.size() < 2) {
    cout << "skipped sequence: " << met.name_ << endl;
    return false;
  }

  // For the reuse of a class instance, clear() is necessary before resize()  
  // because resize() a vector to the same size as the previous use DOES NOT 
  // overwrite its values. 
  // e.g. 
  // EmitProb(NState, vector<ValueType>(100));
  // EmitProb(NState, vector<ValueType>(1000));
  // results in NState x 100 vector rather than NState x 1000 vector
  // because resize() to the same size (i.e. NState) DOES NOT 
  // overwrite their values (i.e. vector<ValueType>(100))
  InitProb.clear();
  TransProb.clear();
  EmitProb.clear();
  vtb_.clear(); 
  trc_vtb_.clear(); 
  path_vtb_.clear(); 
  pdc_.clear(); 
  trc_pdc_.clear(); 
  path_pdc_.clear(); 
  ppr_.clear();
  fwd_.clear();
  bwd_.clear();
  total_ = 0.0;

  // flags
  flg_vtb_ = false;
  flg_trc_vtb_ = false;
  flg_pdc_ = false;
  flg_trc_pdc_ = false;
  flg_ppr_ = false;
  flg_fwd_ = false;
  flg_bwd_ = false;

  // InitProb, TransProb, EmitProb
  reset_param(met, alpha);

  return true;
}

void ProbabilityModel::
print_param(bool logspc) 
{
  cout << "initial distribution:" << endl << flush;
  for (uint i=0; i!=NState; ++i) {
    ValueType val = InitProb[i];
    if (logspc) {
      if (val > NEG_INF) val = Exp(val);
      else val = 0.0;
    }
    cout << val << "\t";
  }
  cout << endl << flush;

  cout << "transition probabilities:" << endl << flush;
  for (uint i=0; i!=NState; ++i) {
    for (uint j=0; j!=NState; ++j) {
      ValueType val = TransProb[i][j];
      if (logspc) {
	if (val > NEG_INF) val = Exp(val);
	else val = 0.0;
      }
      cout << val << "\t";
    }
    cout << endl << flush;
  }
  cout << endl << flush;
}

void ProbabilityModel::
print_table(const vector<vector<ValueType> >& tbl, bool logspc)
{
  cout << "print table:" << endl << flush;
  for (uint i=0; i!=dpsize_; ++i) {
    cout << i << "\t";
    for (uint j=0; j!=NState; ++j) {
      ValueType val = tbl[j][i];
      if (logspc) {
	if (val > NEG_INF) val = Exp(val);
	else val = 0.0;
      }
      cout << val << "\t";
    }
    cout << endl << flush;
  }
  cout << endl << flush;
}

void ProbabilityModel::
viterbi() 
{
  flg_vtb_ = false;
  vtb_.resize(NState, vector<ValueType>(dpsize_, NEG_INF));
  trc_vtb_.resize(NState, vector<uint>(dpsize_, NState));

  for (uint j=0; j!=NState; ++j) 
    vtb_[j][0] = InitProb[j] + EmitProb[j][0];
  for (uint i=1; i!=dpsize_; ++i) {
    for (uint j=0; j!=NState; ++j) {
      ValueType max = NEG_INF;
      uint idx = NState;
      for (uint k=0; k!=NState; ++k) {
	ValueType v = vtb_[k][i-1] + TransProb[k][j] + EmitProb[j][i];
	if (v > max) {
	  max = v;
	  idx = k;
	}
      }
      vtb_[j][i] = max;
      trc_vtb_[j][i] = idx;
    }
  }

  flg_vtb_ = true;
}

void ProbabilityModel::
pdecode(ValueType gcpg, ValueType ggap) 
{
  assert(flg_ppr_);
  flg_pdc_ = false;
  pdc_.resize(NState, vector<ValueType>(dpsize_, NEG_INF));
  trc_pdc_.resize(NState, vector<uint>(dpsize_, NState));

  for (uint j=0; j!=NCPGState; ++j) 
    pdc_[j][0] = (gcpg + 1.0) * ppr_[j][0] - 1.0;

  for (uint i=1; i!=dpsize_; ++i) {
    for (uint j=0; j!=NState; ++j) {
      ValueType max = NEG_INF;
      uint idx = NState;
      for (uint k=0; k!=NState; ++k) {
	if (! (TransProb[k][j] > NEG_INF)) continue;
	ValueType v = pdc_[k][i-1];
	if (j < NCPGState) 
	  v += (gcpg + 1.0) * ppr_[j][i] - 1.0;
	else 
	  v += (ggap + 1.0) * ppr_[j][i] - 1.0; 
	if (v > max) {
	  max = v;
	  idx = k;
	}
      }
      pdc_[j][i] = max;
      trc_pdc_[j][i] = idx;
    }
  }

  flg_pdc_ = true;
}

void ProbabilityModel::
trace(vector<uint>& path, 
      const vector<vector<ValueType> >& tbl,
      const vector<vector<uint> >& trc_tbl) 
{
  path.resize(dpsize_, NState);

  ValueType max = NEG_INF;
  uint idx = NState;
  for (uint j=0; j!=NState; ++j) {
    if (tbl[j][dpsize_-1] > max) {
      max = tbl[j][dpsize_-1];
      idx = j;
    }
  }

  path[dpsize_-1] = idx;
  for (uint i=1; i!=dpsize_; ++i) {
    assert(path[dpsize_-i] != NState);
    path[dpsize_-i-1] = trc_tbl[ path[dpsize_-i] ][dpsize_-i];
  }
}


void ProbabilityModel::
trace_viterbi() 
{
  assert(flg_vtb_);
  flg_trc_vtb_ = false;
  trace(path_vtb_, vtb_, trc_vtb_);
  flg_trc_vtb_ = true;
}

void ProbabilityModel::
trace_pdecode() 
{
  assert(flg_pdc_);
  flg_trc_pdc_ = false;
  trace(path_pdc_, pdc_, trc_pdc_);
  flg_trc_pdc_ = true;
}

void ProbabilityModel::
posterior()
{
  assert(flg_fwd_ && flg_bwd_);
  flg_ppr_ = false;
  ppr_.resize(NState, vector<ValueType>(dpsize_, NEG_INF));

  for (uint i=0; i!=dpsize_; ++i) {
    for (uint j=0; j!=NState; ++j) {
      ValueType val = fwd_[j][i] + bwd_[j][i] - total_;
      if (val > NEG_INF) 
	ppr_[j][i] = Exp(val);
      else 
	ppr_[j][i] = 0.0;
    }
  }

  // round numerical errors
  for (uint i=0; i!=dpsize_; ++i) {
    ValueType sum = 0.0;
    for (uint j=0; j!=NState; ++j) 
      sum += ppr_[j][i];
    for (uint j=0; j!=NState; ++j) 
      ppr_[j][i] /= sum;
  }

  flg_ppr_ = true;
}

void ProbabilityModel::
fwdbwd() 
{
  // forward
  flg_fwd_ = false;
  fwd_.resize(NState, vector<ValueType>(dpsize_, NEG_INF));
  for (uint j=0; j!=NState; ++j) 
    fwd_[j][0] = InitProb[j] + EmitProb[j][0];

  for (uint i=1; i!=dpsize_; ++i) {
    for (uint j=0; j!=NState; ++j) {
      ValueType val = NEG_INF;      
      for (uint k=0; k!=NState; ++k) {
	ValueType v = fwd_[k][i-1] + TransProb[k][j] + EmitProb[j][i];
	if (v > NEG_INF) Fast_LogPlusEquals(val, v);
      }
      fwd_[j][i] = val;
    }
  }
  flg_fwd_ = true;

  // backward
  flg_bwd_ = false;
  bwd_.resize(NState, vector<ValueType>(dpsize_, NEG_INF));
  for (uint j=0; j!=NState; ++j) 
    bwd_[j][dpsize_-1] = Log(1.0);

  for (uint i=1; i!=dpsize_; ++i) {
    for (uint j=0; j!=NState; ++j) {
      ValueType val = NEG_INF;      
      for (uint k=0; k!=NState; ++k) {
	ValueType v = TransProb[j][k] + bwd_[k][dpsize_-i] + EmitProb[k][dpsize_-i];
	if (v > NEG_INF) Fast_LogPlusEquals(val, v);
      }
      bwd_[j][dpsize_-i-1] = val;
    }
  }
  flg_bwd_ = true;

  // total probability
  ValueType total_fwd = NEG_INF;
  ValueType total_bwd = NEG_INF;
  for (uint j=0; j!=NState; ++j) {
    ValueType val = fwd_[j][dpsize_-1] + bwd_[j][dpsize_-1];
    if (val > NEG_INF) Fast_LogPlusEquals(total_fwd, val);
  }
  for (uint j=0; j!=NState; ++j) {
    ValueType val = fwd_[j][0] + bwd_[j][0];
    if (val > NEG_INF) Fast_LogPlusEquals(total_bwd, val);
  }

  total_ = (total_fwd + total_bwd) / 2;
  //cout << "log total probability: " 
  //     << total_ << "\t" << total_fwd << "\t" << total_bwd << endl;
}

void ProbabilityModel::
em(uint nitr, bool verbose) 
{
  const ValueType torelance = 1e-4;

  if (verbose) {
    print_param(true);
    progress("parameter estimation");
  }
  for (uint itr=0; itr!=nitr; ++itr) {
    if (verbose) cout << "round " << itr << endl << flush;
    else cout << "." << flush;

    if (verbose) progress("forward backward");
    fwdbwd();

    if (verbose) progress("expected counts");
    vector<ValueType> ex_init(NState);
    for (uint j=0; j!=NState; ++j) {
      ValueType val = InitProb[j] + EmitProb[j][0] + bwd_[j][0] - total_;
      if (val > NEG_INF) ex_init[j] = Exp(val);
      else ex_init[j] = 0.0;
    }

    vector<vector<ValueType> > ex_trans(NState, vector<ValueType>(NState));
    for (uint j=0; j!=NState; ++j) {
      for (uint k=0; k!=NState; ++k) {
	if (! (TransProb[j][k] > NEG_INF)) { 
	  ex_trans[j][k] = 0.0;
	  continue;
	}
	ValueType val = NEG_INF;
	for (uint i=1; i!=dpsize_; ++i) {
	  ValueType v = fwd_[j][i-1] + TransProb[j][k] + EmitProb[k][i] + bwd_[k][i] - total_;
	  if (v > NEG_INF) Fast_LogPlusEquals(val, v);
	}
	if (val > NEG_INF) ex_trans[j][k] = Exp(val);
	else ex_trans[j][k] = 0.0;
      }
    }

    /*
    cout << "expected counts for initial distribution:" << endl << flush;
    for (uint j=0; j!=NState; ++j) 
      cout << ex_init[j] << "\t";
    cout << endl << flush;
    
    cout << "expected counts for transition probabilities:" << endl << flush;
    for (uint j=0; j!=NState; ++j) {
      for (uint k=0; k!=NState; ++k) 
	cout << ex_trans[j][k] << "\t";
      cout << endl << flush;
    }
    cout << endl << flush;
    */

    if (verbose) progress("parameter update");
    bool conv = true;
    ValueType sum = 0.0;
    for (uint j=0; j!=NState; ++j) 
      sum += ex_init[j];
    for (uint j=0; j!=NState; ++j) {
      if (ex_init[j] > 0.0) {
	ValueType val = ex_init[j] / sum;
	if (fabs(Exp(InitProb[j]) - val) > torelance) conv = false;
	InitProb[j] = Log(val);
      }
      else {
	if (Exp(InitProb[j]) > torelance) conv = false;
	InitProb[j] = NEG_INF;
      }
    }

    for (uint j=0; j!=NState; ++j) {
      sum = 0.0;
      for (uint k=0; k!=NState; ++k) 
	sum += ex_trans[j][k];
      for (uint k=0; k!=NState; ++k) {
	if (ex_trans[j][k] > 0.0) {
	  ValueType val = ex_trans[j][k] / sum;
	  if (fabs(Exp(TransProb[j][k]) - val) > torelance) conv = false;
	  TransProb[j][k] = Log(val);
	}
	else {
	  if (Exp(TransProb[j][k]) > torelance) conv = false;
	  TransProb[j][k] = NEG_INF;
	}
      }
    }
    if (verbose) print_param(true);

    if (conv) {
      progress("convergence");
      break;
    }
  } // for (uint itr=0; itr!=nitr; ++itr) {
  cout << flush;
}

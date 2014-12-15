#include <iostream>
#include <cassert>

#include "SlimProbabilityModel.hh"

using namespace std;

bool SlimProbabilityModel::
reset(const MethylList& met, const GlobalStatistics& gstat, bool noncpg)
{
  // noncpg_
  noncpg_ = noncpg;

  // dpsize_, pathsize_
  dpsize_ = met.pos_size();
  pathsize_ = met.pos_.back() - met.pos_[0] + 1;
  cout << "DP size: " << dpsize_ << endl;
  cout << "path size: " << pathsize_ << endl;
  if (dpsize_ < 2 || pathsize_ < 2) {
    cout << "skipped sequence: " << met.name_ << endl;
    return false;
  }

  // Dist
  Dist.clear();
  Dist.resize(dpsize_-1);
  for (uint i=0; i!=dpsize_-1; ++i) {
    Dist[i] = met.pos_[i+1] - met.pos_[i];
    if (Dist[i] < 2 && (! noncpg_)) {
      cout << "error: distance between neighbor CpGs must not be less than 2," << endl
	   << met.name_ << "\t" << met.pos_[i] << "\t" << met.pos_[i+1] << endl;
      return false;
    }
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
  TransProbDist.clear();
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

  // InitProb, TransProb, EmitProb, TransProbDist,
  reset_param(met, gstat);

  return true;
}

void SlimProbabilityModel::
print_param(bool logspc) 
{
  cout << "initial distribution:" << endl << flush;
  for (uint i=0; i!=NCPGState; ++i) {
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

void SlimProbabilityModel::
print_table(const vector<vector<ValueType> >& tbl, bool logspc)
{
  cout << "print table:" << endl << flush;
  for (uint i=0; i!=dpsize_; ++i) {
    cout << i << "\t";
    for (uint j=0; j!=NCPGState; ++j) {
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

void SlimProbabilityModel::
viterbi() 
{
  flg_vtb_ = false;
  vtb_.resize(NCPGState, vector<ValueType>(dpsize_, NEG_INF));
  trc_vtb_.resize(NCPGState, vector<uint>(dpsize_, NCPGState));

  for (uint j=0; j!=NCPGState; ++j) 
    vtb_[j][0] = InitProb[j] + EmitProb[j][0];
  for (uint i=1; i!=dpsize_; ++i) {
    for (uint j=0; j!=NCPGState; ++j) {
      ValueType max = NEG_INF;
      uint idx = NCPGState;
      uint jdx = NGAPState;
      if (Dist[i-1] == 1) {
	for (uint k=0; k!=NCPGState; ++k) {
	  ValueType v = vtb_[k][i-1] + TransProb[k][j] + EmitProb[j][i]; 
	  if (v > max) {
	    max = v;
	    idx = k;
	  }	  
	}
      }
      else {
	for (uint k=0; k!=NCPGState; ++k) {
	  for (uint l=NCPGState; l!=NState; ++l) {
	    ValueType v = vtb_[k][i-1] + TransProb[k][l] 
	      + TransProbDist[l-NCPGState][i-1] + TransProb[l][j] + EmitProb[j][i]; 
	    if (v > max) {
	      max = v;
	      idx = k;
	      jdx = l;
	    }	  
	  }
	}
      }
      vtb_[j][i] = max;
      trc_vtb_[j][i] = idx * NState + jdx; 
    }
  }

  flg_vtb_ = true;
}

void SlimProbabilityModel::
pdecode(ValueType gcpg, ValueType ggap) 
{
  assert(flg_ppr_);
  flg_pdc_ = false;
  pdc_.resize(NState, vector<ValueType>(dpsize_, NEG_INF));
  trc_pdc_.resize(NState, vector<uint>(dpsize_, NCPGState));

  for (uint j=0; j!=NCPGState; ++j) 
    pdc_[j][0] = (gcpg + 1.0) * ppr_[j][0] - 1.0;
  for (uint i=1; i!=dpsize_; ++i) {
    for (uint j=0; j!=NCPGState; ++j) {
      ValueType max = NEG_INF;
      uint idx = NCPGState;
      uint jdx = NGAPState;
      if (Dist[i-1] == 1) {
	for (uint k=0; k!=NCPGState; ++k) {
	  if (! (TransProb[k][j] > NEG_INF)) continue;
	  ValueType v = pdc_[k][i-1] + (gcpg + 1.0) * ppr_[k][i] - 1.0;; 
	  if (v > max) {
	    max = v;
	    idx = k;
	  }
	}
      }
      else {
	for (uint k=0; k!=NCPGState; ++k) {
	  for (uint l=NCPGState; l!=NState; ++l) {
	    if (! (TransProb[k][l] > NEG_INF && TransProb[l][j] > NEG_INF)) continue;
	    ValueType v = pdc_[k][i-1] 
	      + (ggap + 1.0) * ppr_[l][i-1] - (Dist[i-1] - 1.0)
	      + (gcpg + 1.0) * ppr_[k][i] - 1.0;
	    if (v > max) {
	      max = v;
	      idx = k;
	      jdx = l;
	    }
	  }
	}
      }
      pdc_[j][i] = max;
      trc_vtb_[j][i] = idx * NState + jdx; 
    }
  }

  flg_pdc_ = true;
}

void SlimProbabilityModel::
trace(vector<uint>& path, 
      const vector<vector<ValueType> >& tbl,
      const vector<vector<uint> >& trc_tbl) 
{
  path.resize(pathsize_, NState);

  ValueType max = NEG_INF;
  uint idx = NCPGState;
  uint jdx = NGAPState;
  for (uint j=0; j!=NCPGState; ++j) {
    if (tbl[j][dpsize_-1] > max) {
      max = tbl[j][dpsize_-1];
      idx = j;
    }
  }

  uint p = pathsize_ - 1;
  path[p] = idx;
  for (uint i=1; i!=dpsize_; ++i) {
    uint trc = trc_tbl[ path[p--] ][dpsize_-i];
    idx = trc / NState;
    jdx = trc % NState;
    assert(idx < NCPGState && jdx >= NCPGState && jdx < NState);
    for (uint j=0; j!=Dist[dpsize_-i-1]-1; ++j) 
      path[p--] = jdx;
    path[p] = idx;
  }
}

void SlimProbabilityModel::
trace_viterbi() 
{
  assert(flg_vtb_);
  flg_trc_vtb_ = false;
  trace(path_vtb_, vtb_, trc_vtb_);
  flg_trc_vtb_ = true;
}

void SlimProbabilityModel::
trace_pdecode() 
{
  assert(flg_pdc_);
  flg_trc_pdc_ = false;
  trace(path_pdc_, pdc_, trc_pdc_);
  flg_trc_pdc_ = true;
}

void SlimProbabilityModel::
posterior()
{
  assert(flg_fwd_ && flg_bwd_);
  flg_ppr_ = false;
  ppr_.resize(NState, vector<ValueType>(dpsize_, NEG_INF));

  for (uint i=0; i!=dpsize_; ++i) {
    for (uint j=0; j!=NCPGState; ++j) {
      ValueType val = fwd_[j][i] + bwd_[j][i] - total_;
      if (val > NEG_INF) ppr_[j][i] = Exp(val);
      else ppr_[j][i] = 0.0;
    }
  }
  for (uint i=0; i!=dpsize_-1; ++i) {
    for (uint j=NCPGState; j!=NState; ++j) {
      ValueType val = fwd_[j][i] + TransProbDist[j-NCPGState][i] + bwd_[j][i+1] - total_;
      if (val > NEG_INF) ppr_[j][i] = Exp(val);
      else ppr_[j][i] = 0.0;
    }
  }

  // round numerical errors
  for (uint i=0; i!=dpsize_; ++i) {
    ValueType sum = 0.0;
    for (uint j=0; j!=NCPGState; ++j) 
      sum += ppr_[j][i];
    for (uint j=0; j!=NCPGState; ++j) 
      ppr_[j][i] /= sum;

    sum = 0.0;
    for (uint j=NCPGState; j!=NState; ++j) 
      sum += ppr_[j][i];
    for (uint j=NCPGState; j!=NState; ++j) 
      ppr_[j][i] /= sum;
  }

  flg_ppr_ = true;
}

void SlimProbabilityModel::
fwdbwd() 
{
  // forward
  flg_fwd_ = false;
  fwd_.resize(NState, vector<ValueType>(dpsize_, NEG_INF));
  for (uint j=0; j!=NCPGState; ++j) 
    fwd_[j][0] = InitProb[j] + EmitProb[j][0];

  for (uint i=1; i!=dpsize_; ++i) {
    for (uint j=0; j!=NCPGState; ++j) {
      ValueType val = NEG_INF;
      if (Dist[i-1] == 1) {
	for (uint k=0; k!=NCPGState; ++k) {
	  ValueType v = fwd_[k][i-1] + TransProb[k][j] + EmitProb[j][i];
	  if (v > NEG_INF) Fast_LogPlusEquals(val, v);
	}
      }
      else {
	for (uint k=0; k!=NCPGState; ++k) {
	  for (uint l=NCPGState; l!=NState; ++l) {
	    ValueType v = fwd_[k][i-1] + TransProb[k][l] 
	      + TransProbDist[l-NCPGState][i-1] + TransProb[l][j] + EmitProb[j][i]; 
	    if (v > NEG_INF) Fast_LogPlusEquals(val, v);
	  }
	}
      }
      fwd_[j][i] = val;
    }
    for (uint j=NCPGState; j!=NState; ++j) {
      ValueType val = NEG_INF;
      for (uint k=0; k!=NCPGState; ++k) {
	ValueType v = fwd_[k][i-1] + TransProb[k][j]; 
	if (v > NEG_INF) Fast_LogPlusEquals(val, v);
      }
      fwd_[j][i-1] = val;
    }
  }
  flg_fwd_ = true;
  //print_table(fwd_, false);

  // backward
  flg_bwd_ = false;
  bwd_.resize(NState, vector<ValueType>(dpsize_, NEG_INF));
  for (uint j=0; j!=NCPGState; ++j) 
    bwd_[j][dpsize_-1] = Log(1.0);

  for (uint i=1; i!=dpsize_; ++i) {
    for (uint j=0; j!=NCPGState; ++j) {
      ValueType val = NEG_INF;      
      if (Dist[dpsize_-i-1] == 1) {
	for (uint k=0; k!=NCPGState; ++k) {
	  ValueType v = TransProb[j][k] + bwd_[k][dpsize_-i] + EmitProb[k][dpsize_-i];
	  if (v > NEG_INF) Fast_LogPlusEquals(val, v);
	}
      }
      else {
	for (uint k=0; k!=NCPGState; ++k) {
	  for (uint l=NCPGState; l!=NState; ++l) {
	    ValueType v = TransProb[j][l] + TransProbDist[l-NCPGState][dpsize_-i-1] 
	      + TransProb[l][k] + bwd_[k][dpsize_-i] + EmitProb[k][dpsize_-i];
	    if (v > NEG_INF) Fast_LogPlusEquals(val, v);
	  }
	}
      }
      bwd_[j][dpsize_-i-1] = val;
    }
    for (uint j=NCPGState; j!=NState; ++j) {
      ValueType val = NEG_INF;      
      for (uint k=0; k!=NCPGState; ++k) {
	ValueType v = TransProb[j][k] + bwd_[k][dpsize_-i] + EmitProb[k][dpsize_-i];
	if (v > NEG_INF) Fast_LogPlusEquals(val, v);
      }
      bwd_[j][dpsize_-i] = val;
    }
  }
  flg_bwd_ = true;

  // total probability
  ValueType total_fwd = NEG_INF;
  ValueType total_bwd = NEG_INF;
  for (uint j=0; j!=NCPGState; ++j) {
    ValueType val = fwd_[j][dpsize_-1] + bwd_[j][dpsize_-1];
    if (val > NEG_INF) Fast_LogPlusEquals(total_fwd, val);
  }
  for (uint j=0; j!=NCPGState; ++j) {
    ValueType val = fwd_[j][0] + bwd_[j][0];
    if (val > NEG_INF) Fast_LogPlusEquals(total_bwd, val);
  }

  total_ = (total_fwd + total_bwd) / 2;
  //cout << "log total probability: " 
  //     << total_ << "\t" << total_fwd << "\t" << total_bwd << endl;
}

void SlimProbabilityModel::
em(uint nitr, bool verbose) 
{
  const ValueType tolerance = 1e-4;

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
    vector<ValueType> ex_init(NCPGState);
    for (uint j=0; j!=NCPGState; ++j) {
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
	if (j < NCPGState && k < NCPGState) { // CpG to CpG
	  for (uint i=0; i!=dpsize_-1; ++i) {
	    if (Dist[i] == 1) {
	      ValueType v = fwd_[j][i] + TransProb[j][k] + bwd_[k][i+1] + EmitProb[k][i+1] - total_;
	      if (v > NEG_INF) Fast_LogPlusEquals(val, v);
	    }
	  }
	}
	else if (j < NCPGState && k >= NCPGState) { // CpG to GAP
	  for (uint i=0; i!=dpsize_-1; ++i) {
	    ValueType v = fwd_[j][i] + TransProb[j][k] 
	      + TransProbDist[k-NCPGState][i] + bwd_[k][i+1] - total_;
	    if (v > NEG_INF) Fast_LogPlusEquals(val, v);
	  }
	}
	else if (j >= NCPGState && k < NCPGState) { // GAP to CpG
	  for (uint i=0; i!=dpsize_-1; ++i) {
	    ValueType v = fwd_[j][i] + TransProbDist[j-NCPGState][i]
	      + TransProb[j][k] + bwd_[k][i+1] + EmitProb[k][i+1] - total_;
	    if (v > NEG_INF) Fast_LogPlusEquals(val, v);
	  }
	}
	else if (j==k) { // GAP loop
	  for (uint i=0; i!=dpsize_-1; ++i) {
	    ValueType v = fwd_[j][i] + TransProbDist[j-NCPGState][i] + bwd_[j][i+1] - total_;
	    v += Log((ValueType) Dist[i] - 2);
	    if (v > NEG_INF) Fast_LogPlusEquals(val, v);
	  }
	}
	else {
	  assert(true);
	}
	if (val > NEG_INF) ex_trans[j][k] = Exp(val);
	else ex_trans[j][k] = 0.0;
      }
    }

    /*
    cout << "expected counts for initial distribution:" << endl << flush;
    for (uint j=0; j!=NCPGState; ++j) 
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
    for (uint j=0; j!=NCPGState; ++j) 
      sum += ex_init[j];
    for (uint j=0; j!=NCPGState; ++j) {
      if (ex_init[j] > 0.0) {
	ValueType val = ex_init[j] / sum;
	if (fabs(Exp(InitProb[j]) - val) > tolerance) conv = false;
	InitProb[j] = Log(val);
      }
      else {
	if (Exp(InitProb[j]) > tolerance) conv = false;
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
	  if (fabs(Exp(TransProb[j][k]) - val) > tolerance) conv = false;
	  TransProb[j][k] = Log(val);
	}
	else {
	  if (Exp(TransProb[j][k]) > tolerance) conv = false;
	  TransProb[j][k] = NEG_INF;
	}
      }
    }
    if (verbose) print_param(true);

    for (uint i=0; i!=dpsize_-1; ++i) {
      for (uint j=NCPGState; j!=NState; ++j) {
	if (Dist[i] == 1) TransProbDist[j-NCPGState][i] = NEG_INF;
	else TransProbDist[j-NCPGState][i] = (Dist[i] - 2) * TransProb[j][j];
      }
    }

    if (conv) {
      progress("convergence");
      break;
    }
  } // for (uint itr=0; itr!=nitr; ++itr) {
  cout << flush;
}

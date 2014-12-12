#include <algorithm>
#include <iostream>
#include <cstdlib>

#include "GlobalStatistics.hh"

using namespace std;

bool GlobalStatistics::
reset(const vector<MethylList>& data, 
      uint nmix, bool nobeta, uint nitr, uint nsmp, bool verbose)
{
  const uint precomp_size = 2000;

  gamma_size_ = precomp_size;
  gamma_.clear();
  gamma_.resize(gamma_size_);
  gamma_[0] = 0.0;
  for (uint i=1; i<gamma_size_; ++i) 
    gamma_[i] = fast_gamma((double) i);

  digamma_size_ = precomp_size;
  digamma_.clear();
  digamma_.resize(digamma_size_);
  digamma_[0] = 0.0;
  for (uint i=1; i<digamma_size_; ++i) 
    digamma_[i] = fast_digamma((double) i);

  nmix_ = nmix;
  nobeta_ = nobeta;
  if (nobeta_) return true;
  reset_beta_mixture(data, 1, nitr, nsmp, verbose);
  reset_beta_mixture(data, 2, nitr, nsmp, verbose);
  if (nmix_ == 0) nmix_ = 1;

  print_param(1);
  print_param(2);

  return true;
}

GlobalStatistics::ValueType GlobalStatistics::
pup(const MethylList& met, uint i) const
{
  if (nobeta_) return pup_pseudo_count(met, i);
  else return pup_beta_mixture(met, i);
}

GlobalStatistics::ValueType GlobalStatistics::
pdown(const MethylList& met, uint i) const
{
  if (nobeta_) return pdown_pseudo_count(met, i);
  else return pdown_beta_mixture(met, i);
}

GlobalStatistics::ValueType GlobalStatistics::
pnochange(const MethylList& met, uint i) const
{
  if (nobeta_) return pnochange_pseudo_count(met, i);
  else return pnochange_beta_mixture(met, i);
}

GlobalStatistics::ValueType GlobalStatistics::
pup(float m1, float u1, float m2, float u2) const
{
  if (nobeta_) return pup_pseudo_count(m1, u1, m2, u2);
  else return pup_beta_mixture(m1, u1, m2, u2);
}

GlobalStatistics::ValueType GlobalStatistics::
pdown(float m1, float u1, float m2, float u2) const
{
  if (nobeta_) return pdown_pseudo_count(m1, u1, m2, u2);
  else return pdown_beta_mixture(m1, u1, m2, u2);
}

GlobalStatistics::ValueType GlobalStatistics::
pnochange(float m1, float u1, float m2, float u2) const
{
  if (nobeta_) return pnochange_pseudo_count(m1, u1, m2, u2);
  else return pnochange_beta_mixture(m1, u1, m2, u2);
}

void GlobalStatistics::
print_param(uint c)
{
  assert(c == 1 || c == 2);

  vector<ValueType>* whp;
  vector<float>* ahp;
  vector<float>* bhp;
  vector<ValueType>* wlp;
  vector<float>* alp;
  vector<float>* blp;
  if (c == 1) {
    whp = &wh1_;
    ahp = &ah1_;
    bhp = &bh1_;
    wlp = &wl1_;
    alp = &al1_;
    blp = &bl1_;
  }
  else {
    whp = &wh2_;
    ahp = &ah2_;
    bhp = &bh2_;
    wlp = &wl2_;
    alp = &al2_;
    blp = &bl2_;
  }

  cout << "beta mixture for sample " << c << ":" << endl << flush;
  for (uint i=0; i!=(*whp).size(); ++i) 
    cout << "high\t" << (*whp)[i] << "\t" << (*ahp)[i] << "\t" << (*bhp)[i] << endl << flush;
  for (uint i=0; i!=(*wlp).size(); ++i) 
    cout << "low\t"  << (*wlp)[i] << "\t" << (*alp)[i] << "\t" << (*blp)[i] << endl << flush;
  cout << endl << flush;
}

struct IJR 
{
  uint i_, j_, r_;

  IJR(uint i, uint j, uint r) : i_(i), j_(j), r_(r) {}
  IJR(const IJR& x) : i_(x.i_), j_(x.j_), r_(x.r_) {}
  ~IJR() {}
};

void GlobalStatistics::
reset_beta_mixture(const vector<MethylList>& data, uint c, 
		   uint nitr, uint nsmp, bool verbose)
{
  const double tolerance = 1e-2;
  const double regcoef = 1e-4;
  //const double regcoef = 1e-6;
  const uint nitr_outer = nitr;
  const uint nitr_inner = nitr;
  assert(c == 1 || c == 2);

  vector<ValueType>* whp;
  vector<float>* ahp;
  vector<float>* bhp;
  vector<ValueType>* wlp;
  vector<float>* alp;
  vector<float>* blp;
  const float* mp = NULL;
  const float* up = NULL;
  const float* np = NULL;
  if (c == 1) {
    whp = &wh1_;
    ahp = &ah1_;
    bhp = &bh1_;
    wlp = &wl1_;
    alp = &al1_;
    blp = &bl1_;
  }
  else {
    whp = &wh2_;
    ahp = &ah2_;
    bhp = &bh2_;
    wlp = &wl2_;
    alp = &al2_;
    blp = &bl2_;
  }

  (*whp).clear();
  (*ahp).clear();
  (*bhp).clear();
  (*wlp).clear();
  (*alp).clear();
  (*blp).clear();

  if (nmix_ == 0) {
    (*whp).resize(2);
    (*ahp).resize(2);
    (*bhp).resize(2);
    (*wlp).resize(2);
    (*alp).resize(2);
    (*blp).resize(2);
    // fixed one-sided beta
    (*whp)[0] = 0.1;
    (*ahp)[0] = 8.0;
    (*bhp)[0] = 1.0;
    (*wlp)[0] = 0.1;
    (*alp)[0] = 1.0;
    (*blp)[0] = 8.0;
    // uniform
    (*whp)[1] = 0.4;
    (*ahp)[1] = 1.0;
    (*bhp)[1] = 1.0;
    (*wlp)[1] = 0.4;
    (*alp)[1] = 1.0;
    (*blp)[1] = 1.0;
    return;
  }

  srand(RAND_SEED);
  (*whp).resize(nmix_+1);
  (*ahp).resize(nmix_+1);
  (*bhp).resize(nmix_+1);
  (*wlp).resize(nmix_+1);
  (*alp).resize(nmix_+1);
  (*blp).resize(nmix_+1);
  for (uint k=0; k!=nmix_; ++k) {
    (*whp)[k] = 0.1 / nmix_;
    (*ahp)[k] = 1.0 + (10.0 * rand() / RAND_MAX);
    (*bhp)[k] = 1.0;
    (*wlp)[k] = 0.1 / nmix_;
    (*alp)[k] = 1.0; 
    (*blp)[k] = 1.0 + (10.0 * rand() / RAND_MAX);
  }
  (*whp)[nmix_] = 0.4;
  (*ahp)[nmix_] = 1.0;
  (*bhp)[nmix_] = 1.0;
  (*wlp)[nmix_] = 0.4;
  (*alp)[nmix_] = 1.0; 
  (*blp)[nmix_] = 1.0;

  srand(RAND_SEED);
  vector<IJR> ijrdx;
  uint npos = 0;
  double scale = 0.0;
  for (uint i=0; i!=data.size(); ++i) {
    for (uint j=0; j!=data[i].pos_size(); ++j) {
      for (uint r=0; r!=data[i].rep_size(c, j); ++r) {
	ijrdx.push_back(IJR(i,j,r));
	npos++;
      }
    }
  }
  if (npos > nsmp) {
    random_shuffle(ijrdx.begin(), ijrdx.end());
    npos = nsmp;
  }
  scale = 1.0 / npos;

  // posterior probabilities p(component k | data[i][j], w_, a_, b_) 
  // likelihood p(data[i][j] | component k, w_, a_, b_) 
  vector<vector<double> > pstrh(npos, vector<double>(nmix_+1, 0.0));
  vector<vector<double> > lkhdh(npos, vector<double>(nmix_+1, 0.0));
  vector<vector<double> > pstrl(npos, vector<double>(nmix_+1, 0.0));
  vector<vector<double> > lkhdl(npos, vector<double>(nmix_+1, 0.0));

  // precomputation for fast_gamma
  vector<double> gm_m(npos);
  vector<double> gm_u(npos);
  vector<double> gm_n(npos);
  for (uint ijr=0; ijr!=npos; ++ijr) {
    uint i = ijrdx[ijr].i_;
    uint j = ijrdx[ijr].j_;
    uint r = ijrdx[ijr].r_;
    data[i].get_pointer(mp, up, np, c, j, r);
    gm_m[ijr] = gamma((*mp) + 1.0F);
    gm_u[ijr] = gamma((*up) + 1.0F);
    gm_n[ijr] = gamma((*np) + 1.0F);
  }

  vector<double> ln_ah(nmix_);
  vector<double> ln_bl(nmix_); 
  for (uint k=0; k!=nmix_; ++k) {
    ln_ah[k] = Log((*ahp)[k] - 1.0);
    ln_bl[k] = Log((*blp)[k] - 1.0);
  }

  if (verbose) {
    print_param(c);
    progress("prior estimation");
  }
  for (uint itr_outer=0; itr_outer!=nitr_outer; ++itr_outer) {
    bool conv = true;

    for (uint itr_inner=0; itr_inner!=nitr_inner; ++itr_inner) {
      if (verbose) cout << "outer round " << itr_outer << " inner ab round " << itr_inner << endl << flush;
      else cout << "." << flush;
      
      vector<float> abh(nmix_+1);
      vector<double> diabh(nmix_+1);
      vector<double> gmabh(nmix_+1);
      vector<double> sum_ah(nmix_+1);
      vector<double> sum_ph(nmix_+1);
      vector<float> abl(nmix_+1);
      vector<double> diabl(nmix_+1);
      vector<double> gmabl(nmix_+1);
      vector<double> sum_bl(nmix_+1);
      vector<double> sum_pl(nmix_+1);
      for (uint k=0; k!=nmix_+1; ++k) {
	abh[k] = (*ahp)[k] + (*bhp)[k];
	diabh[k] = digamma(abh[k]);
	gmabh[k] = gamma(abh[k]);
	sum_ah[k] = 0.0;
	sum_ph[k] = 0.0;
	abl[k] = (*alp)[k] + (*blp)[k];
	diabl[k] = digamma(abl[k]);
	gmabl[k] = gamma(abl[k]);
	sum_bl[k] = 0.0;
	sum_pl[k] = 0.0;
      }
      for (uint ijr=0; ijr!=npos; ++ijr) {
	uint i = ijrdx[ijr].i_;
	uint j = ijrdx[ijr].j_;
	uint r = ijrdx[ijr].r_;
	data[i].get_pointer(mp, up, np, c, j, r);
	vector<double> diabhn(nmix_+1);
	vector<double> nmrth(nmix_+1);
	vector<double> diabln(nmix_+1);
	vector<double> nmrtl(nmix_+1);
	double dnmt = 0.0;
	for (uint k=0; k!=nmix_+1; ++k) {
	  diabhn[k] = digamma(abh[k] + (*np));
	  diabln[k] = digamma(abl[k] + (*np));
	}
	for (uint k=0; k!=nmix_+1; ++k) {
	  double vh = gmabh[k] + gm_n[ijr] - gamma(abh[k] + (*np))
	    + gamma((*ahp)[k] + (*mp)) - gamma((*ahp)[k]) - gm_m[ijr] 
	    + gamma((*bhp)[k] + (*up)) - gamma((*bhp)[k]) - gm_u[ijr];
	  lkhdh[ijr][k] = Exp(vh);
	  nmrth[k] = (*whp)[k] * lkhdh[ijr][k];
	  double vl = gmabl[k] + gm_n[ijr] - gamma(abl[k] + (*np))
	    + gamma((*alp)[k] + (*mp)) - gamma((*alp)[k]) - gm_m[ijr] 
	    + gamma((*blp)[k] + (*up)) - gamma((*blp)[k]) - gm_u[ijr];
	  lkhdl[ijr][k] = Exp(vl);
	  nmrtl[k] = (*wlp)[k] * lkhdl[ijr][k];
	  dnmt += nmrth[k];
	  dnmt += nmrtl[k];
	}
	for (uint k=0; k!=nmix_+1; ++k) {
	  pstrh[ijr][k] = nmrth[k] / dnmt;
	  pstrl[ijr][k] = nmrtl[k] / dnmt;
	}
	for (uint k=0; k!=nmix_+1; ++k) {
	  sum_ah[k] += pstrh[ijr][k] * (digamma((*ahp)[k] + (*mp)) - diabhn[k]);
	  sum_ph[k] += pstrh[ijr][k];
	  sum_bl[k] += pstrl[ijr][k] * (digamma((*blp)[k] + (*up)) - diabln[k]);
	  sum_pl[k] += pstrl[ijr][k];
	}
      }
      
      if (verbose) progress("prior ab update");
      bool conv_ab = true;
      for (uint k=0; k!=nmix_; ++k) {
	// L2 regularizer for preventing large a, b
	ln_ah[k] += scale * (*ahp)[k] * (sum_ph[k] * (diabh[k] - digamma((*ahp)[k])) + sum_ah[k])
	  - regcoef * (*ahp)[k] * ((*ahp)[k] - 1.0);
	ln_bl[k] += scale * (*blp)[k] * (sum_pl[k] * (diabl[k] - digamma((*blp)[k])) + sum_bl[k])
	  - regcoef * (*blp)[k] * ((*blp)[k] - 1.0);
	double val_ah = Exp(ln_ah[k]) + 1.0;
	double val_bl = Exp(ln_bl[k]) + 1.0;
	if (fabs((*ahp)[k] - val_ah) > tolerance) conv_ab = false;
	if (fabs((*blp)[k] - val_bl) > tolerance) conv_ab = false;
	(*ahp)[k] = val_ah;
	(*blp)[k] = val_bl;
      }
      if (verbose) print_param(c);
      
      if (itr_inner > 0) conv = false;
      if (conv_ab) {
	if (itr_inner > 0) conv = false;
	progress("inner ab convergence");
	break;
      }
    } // for (uint itr_inner=0; itr_inner!=nitr_inner; ++itr_inner) {

    for (uint itr_inner=0; itr_inner!=nitr_inner; ++itr_inner) {
      if (verbose) cout << "outer round " << itr_outer << " inner w round " << itr_inner << endl << flush;
      else cout << "." << flush;
      
      vector<double> sum_ph(nmix_+1, 0.0);
      vector<double> sum_pl(nmix_+1, 0.0);
      for (uint ijr=0; ijr!=npos; ++ijr) {
	vector<double> nmrth(nmix_+1);
	vector<double> nmrtl(nmix_+1);
	double dnmt = 0.0;
	for (uint k=0; k!=nmix_+1; ++k) {
	  nmrth[k] = (*whp)[k] * lkhdh[ijr][k];
	  nmrtl[k] = (*wlp)[k] * lkhdl[ijr][k];
	  dnmt += nmrth[k];
	  dnmt += nmrtl[k];
	}
	for (uint k=0; k!=nmix_+1; ++k) {
	  pstrh[ijr][k] = nmrth[k] / dnmt;
	  pstrl[ijr][k] = nmrtl[k] / dnmt;
	}
	for (uint k=0; k!=nmix_+1; ++k) {
	  sum_ph[k] += pstrh[ijr][k];
	  sum_pl[k] += pstrl[ijr][k];
	}
      }

      if (verbose) progress("prior ab update");
      bool conv_w = true;
      double norm = 0.0;
      for (uint k=0; k!=nmix_+1; ++k) {
	norm += sum_ph[k];
	norm += sum_pl[k];
      }
      for (uint k=0; k!=nmix_+1; ++k) {	
	double val_wh = sum_ph[k] / norm;
	double val_wl = sum_pl[k] / norm;
	if (fabs((*whp)[k] - val_wh) > tolerance) conv_w = false;
	if (fabs((*wlp)[k] - val_wl) > tolerance) conv_w = false;
	(*whp)[k] = val_wh;
	(*wlp)[k] = val_wl;
      }
      if (verbose) print_param(c);
      
      if (itr_inner > 0) conv = false;
      if (conv_w) {
	progress("inner w convergence");
	break;
      }
    } // for (uint itr_inner=0; itr_inner!=nitr_inner; ++itr_inner) {
    
    if (conv) {
      progress("outer convergence");
      break;
    }
  } // for (uint itr_outer=0; itr_outer!=nitr_outer; ++itr_outer) {
}

GlobalStatistics::ValueType GlobalStatistics::
pup_beta_mixture(const MethylList& met, uint i) const
{
  const std::vector<float>* m1p = NULL;
  const std::vector<float>* u1p = NULL; 
  const std::vector<float>* n1p = NULL;
  const float* m1sump = NULL;
  const float* u1sump = NULL; 
  const float* n1sump = NULL; 
  const std::vector<float>* m2p = NULL;
  const std::vector<float>* u2p = NULL; 
  const std::vector<float>* n2p = NULL;
  const float* m2sump = NULL;
  const float* u2sump = NULL;
  const float* n2sump = NULL;

  met.get_pointer(m1p, u1p, n1p, m1sump, u1sump, n1sump, 1, i);
  met.get_pointer(m2p, u2p, n2p, m2sump, u2sump, n2sump, 2, i);

  return pup_beta_mixture(*m1sump, *u1sump, *m2sump, *u2sump);
}

GlobalStatistics::ValueType GlobalStatistics::
pdown_beta_mixture(const MethylList& met, uint i) const
{
  const std::vector<float>* m1p = NULL;
  const std::vector<float>* u1p = NULL; 
  const std::vector<float>* n1p = NULL;
  const float* m1sump = NULL;
  const float* u1sump = NULL; 
  const float* n1sump = NULL; 
  const std::vector<float>* m2p = NULL;
  const std::vector<float>* u2p = NULL; 
  const std::vector<float>* n2p = NULL;
  const float* m2sump = NULL;
  const float* u2sump = NULL;
  const float* n2sump = NULL;

  met.get_pointer(m1p, u1p, n1p, m1sump, u1sump, n1sump, 1, i);
  met.get_pointer(m2p, u2p, n2p, m2sump, u2sump, n2sump, 2, i);

  return pdown_beta_mixture(*m1sump, *u1sump, *m2sump, *u2sump);
}

GlobalStatistics::ValueType GlobalStatistics::
pnochange_beta_mixture(const MethylList& met, uint i) const
{
  const std::vector<float>* m1p = NULL;
  const std::vector<float>* u1p = NULL; 
  const std::vector<float>* n1p = NULL;
  const float* m1sump = NULL;
  const float* u1sump = NULL; 
  const float* n1sump = NULL; 
  const std::vector<float>* m2p = NULL;
  const std::vector<float>* u2p = NULL; 
  const std::vector<float>* n2p = NULL;
  const float* m2sump = NULL;
  const float* u2sump = NULL;
  const float* n2sump = NULL;

  met.get_pointer(m1p, u1p, n1p, m1sump, u1sump, n1sump, 1, i);
  met.get_pointer(m2p, u2p, n2p, m2sump, u2sump, n2sump, 2, i);

  return pnochange_beta_mixture(*m1sump, *u1sump, *m2sump, *u2sump);
}

GlobalStatistics::ValueType GlobalStatistics::
pup_pseudo_count(const MethylList& met, uint i) const
{
  const std::vector<float>* m1p = NULL;
  const std::vector<float>* u1p = NULL; 
  const std::vector<float>* n1p = NULL;
  const float* m1sump = NULL;
  const float* u1sump = NULL; 
  const float* n1sump = NULL; 
  const std::vector<float>* m2p = NULL;
  const std::vector<float>* u2p = NULL; 
  const std::vector<float>* n2p = NULL;
  const float* m2sump = NULL;
  const float* u2sump = NULL;
  const float* n2sump = NULL;

  met.get_pointer(m1p, u1p, n1p, m1sump, u1sump, n1sump, 1, i);
  met.get_pointer(m2p, u2p, n2p, m2sump, u2sump, n2sump, 2, i);

  return pup_pseudo_count(*m1sump, *u1sump, *m2sump, *u2sump);
}

GlobalStatistics::ValueType GlobalStatistics::
pdown_pseudo_count(const MethylList& met, uint i) const
{
  const std::vector<float>* m1p = NULL;
  const std::vector<float>* u1p = NULL; 
  const std::vector<float>* n1p = NULL;
  const float* m1sump = NULL;
  const float* u1sump = NULL; 
  const float* n1sump = NULL; 
  const std::vector<float>* m2p = NULL;
  const std::vector<float>* u2p = NULL; 
  const std::vector<float>* n2p = NULL;
  const float* m2sump = NULL;
  const float* u2sump = NULL;
  const float* n2sump = NULL;

  met.get_pointer(m1p, u1p, n1p, m1sump, u1sump, n1sump, 1, i);
  met.get_pointer(m2p, u2p, n2p, m2sump, u2sump, n2sump, 2, i);

  return pdown_pseudo_count(*m1sump, *u1sump, *m2sump, *u2sump);
}

GlobalStatistics::ValueType GlobalStatistics::
pnochange_pseudo_count(const MethylList& met, uint i) const
{
  const std::vector<float>* m1p = NULL;
  const std::vector<float>* u1p = NULL; 
  const std::vector<float>* n1p = NULL;
  const float* m1sump = NULL;
  const float* u1sump = NULL; 
  const float* n1sump = NULL; 
  const std::vector<float>* m2p = NULL;
  const std::vector<float>* u2p = NULL; 
  const std::vector<float>* n2p = NULL;
  const float* m2sump = NULL;
  const float* u2sump = NULL;
  const float* n2sump = NULL;

  met.get_pointer(m1p, u1p, n1p, m1sump, u1sump, n1sump, 1, i);
  met.get_pointer(m2p, u2p, n2p, m2sump, u2sump, n2sump, 2, i);

  return pnochange_pseudo_count(*m1sump, *u1sump, *m2sump, *u2sump);
}


GlobalStatistics::ValueType GlobalStatistics::
pup_beta_mixture(float m1, float u1, float m2, float u2) const
{
  const bool posterior = false;

  const vector<ValueType>* w1p = &(wh1_);
  const vector<float>* a1p = &(ah1_);
  const vector<float>* b1p = &(bh1_);
  const vector<ValueType>* w2p = &(wl2_);
  const vector<float>* a2p = &(al2_);
  const vector<float>* b2p = &(bl2_);

  vector<ValueType> pstr1(nmix_);
  vector<ValueType> pstr2(nmix_);
  vector<ValueType> nmrt;
  ValueType dnmt;

  nmrt.clear();
  nmrt.resize(nmix_);
  dnmt = 0.0;
  for (uint i=0; i!=nmix_; ++i) {
    if (posterior) {
      ValueType v = gamma((*a1p)[i] + (*b1p)[i]) + gamma(m1 + u1 + 1.0F) - gamma((*a1p)[i] + (*b1p)[i] + m1 + u1)
	+ gamma((*a1p)[i] + m1) - gamma((*a1p)[i]) - gamma(m1 + 1.0F) 
	+ gamma((*b1p)[i] + u1) - gamma((*b1p)[i]) - gamma(u1 + 1.0F);
      nmrt[i] = (*w1p)[i] * Exp(v);
    }
    else {
      nmrt[i] = (*w1p)[i];
    }
    dnmt += nmrt[i];
  }
  for (uint i=0; i!=nmix_; ++i) 
    pstr1[i] = nmrt[i] / dnmt;

  nmrt.clear();
  nmrt.resize(nmix_);
  dnmt = 0.0;
  for (uint i=0; i!=nmix_; ++i) {
    if (posterior) {
      ValueType v = gamma((*a2p)[i] + (*b2p)[i]) + gamma(m2 + u2 + 1.0F) - gamma((*a2p)[i] + (*b2p)[i] + m2 + u2)
	+ gamma((*a2p)[i] + m2) - gamma((*a2p)[i]) - gamma(m2 + 1.0F) 
	+ gamma((*b2p)[i] + u2) - gamma((*b2p)[i]) - gamma(u2 + 1.0F);
      nmrt[i] = (*w2p)[i] * Exp(v);
    }
    else {
      nmrt[i] = (*w2p)[i];
    }
    dnmt += nmrt[i];
  }
  for (uint i=0; i!=nmix_; ++i) 
    pstr2[i] = nmrt[i] / dnmt;

  ValueType val1 = NEG_INF;
  ValueType val2 = NEG_INF;
  for (uint i=0; i!=nmix_; ++i) {
    ValueType v;
    if (posterior) v = Log(pstr1[i]) 
      + beta(m1+m1+(*a1p)[i], u1+u1+(*b1p)[i]) - beta(m1+(*a1p)[i], u1+(*b1p)[i]);
    else v = Log(pstr1[i]) 
      + beta(m1+(*a1p)[i], u1+(*b1p)[i]) - beta((*a1p)[i], (*b1p)[i]);
    Fast_LogPlusEquals(val1, v);
  }
  for (uint i=0; i!=nmix_; ++i) {
    ValueType v;
    if (posterior) v = Log(pstr2[i]) 
      + beta(m2+m2+(*a2p)[i], u2+u2+(*b2p)[i]) - beta(m2+(*a2p)[i], u2+(*b2p)[i]);
    else v = Log(pstr2[i]) 
      + beta(m2+(*a2p)[i], u2+(*b2p)[i]) - beta((*a2p)[i], (*b2p)[i]);
    Fast_LogPlusEquals(val2, v);
  }

  ValueType p = comb(m1+u1, m1) + val1 + comb(m2+u2, m2) + val2;
  return p;
}

GlobalStatistics::ValueType GlobalStatistics::
pdown_beta_mixture(float m1, float u1, float m2, float u2) const
{
  const bool posterior = false;

  const vector<ValueType>* w1p = &(wl1_);
  const vector<float>* a1p = &(al1_);
  const vector<float>* b1p = &(bl1_);
  const vector<ValueType>* w2p = &(wh2_);
  const vector<float>* a2p = &(ah2_);
  const vector<float>* b2p = &(bh2_);

  vector<ValueType> pstr1(nmix_);
  vector<ValueType> pstr2(nmix_);
  vector<ValueType> nmrt;
  ValueType dnmt;

  nmrt.clear();
  nmrt.resize(nmix_);
  dnmt = 0.0;
  for (uint i=0; i!=nmix_; ++i) {
    if (posterior) {
      ValueType v = gamma((*a1p)[i] + (*b1p)[i]) + gamma(m1 + u1 + 1.0F) - gamma((*a1p)[i] + (*b1p)[i] + m1 + u1)
	+ gamma((*a1p)[i] + m1) - gamma((*a1p)[i]) - gamma(m1 + 1.0F) 
	+ gamma((*b1p)[i] + u1) - gamma((*b1p)[i]) - gamma(u1 + 1.0F);
      nmrt[i] = (*w1p)[i] * Exp(v);
    }
    else {
      nmrt[i] = (*w1p)[i];
    }
    dnmt += nmrt[i];
  }
  for (uint i=0; i!=nmix_; ++i) 
    pstr1[i] = nmrt[i] / dnmt;

  nmrt.clear();
  nmrt.resize(nmix_);
  dnmt = 0.0;
  for (uint i=0; i!=nmix_; ++i) {
    if (posterior) {
      ValueType v = gamma((*a2p)[i] + (*b2p)[i]) + gamma(m2 + u2 + 1.0F) - gamma((*a2p)[i] + (*b2p)[i] + m2 + u2)
	+ gamma((*a2p)[i] + m2) - gamma((*a2p)[i]) - gamma(m2 + 1.0F) 
	+ gamma((*b2p)[i] + u2) - gamma((*b2p)[i]) - gamma(u2 + 1.0F);
      nmrt[i] = (*w2p)[i] * Exp(v);
    }
    else {
      nmrt[i] = (*w2p)[i];
    }
    dnmt += nmrt[i];
  }
  for (uint i=0; i!=nmix_; ++i) 
    pstr2[i] = nmrt[i] / dnmt;

  ValueType val1 = NEG_INF;
  ValueType val2 = NEG_INF;
  for (uint i=0; i!=nmix_; ++i) {
    ValueType v;
    if (posterior) v = Log(pstr1[i]) 
      + beta(m1+m1+(*a1p)[i], u1+u1+(*b1p)[i]) - beta(m1+(*a1p)[i], u1+(*b1p)[i]);
    else v = Log(pstr1[i]) 
      + beta(m1+(*a1p)[i], u1+(*b1p)[i]) - beta((*a1p)[i], (*b1p)[i]);
    Fast_LogPlusEquals(val1, v);
  }
  for (uint i=0; i!=nmix_; ++i) {
    ValueType v;
    if (posterior) v = Log(pstr2[i]) 
      + beta(m2+m2+(*a2p)[i], u2+u2+(*b2p)[i]) - beta(m2+(*a2p)[i], u2+(*b2p)[i]);
    else v = Log(pstr2[i]) 
      + beta(m2+(*a2p)[i], u2+(*b2p)[i]) - beta((*a2p)[i], (*b2p)[i]);
    Fast_LogPlusEquals(val2, v);
  }

  ValueType p = comb(m1+u1, m1) + val1 + comb(m2+u2, m2) + val2;
  return p;
}

GlobalStatistics::ValueType GlobalStatistics::
pnochange_beta_mixture(float m1, float u1, float m2, float u2) const
{
  const bool posterior = false;

  assert(wh1_.size() == wl1_.size());
  assert(wh2_.size() == wl2_.size());

  vector<ValueType> pstrh1(wh1_.size());
  vector<ValueType> pstrl1(wh1_.size());
  vector<ValueType> pstrh2(wh2_.size());
  vector<ValueType> pstrl2(wh2_.size());
  vector<ValueType> nmrth;
  vector<ValueType> nmrtl;
  ValueType dnmt;

  nmrth.clear();
  nmrth.resize(wh1_.size());
  nmrtl.clear();
  nmrtl.resize(wh1_.size());
  dnmt = 0.0;
  for (uint i=0; i!=wh1_.size(); ++i) {
    if (posterior) {
      /*
      cout << "m1:u1:m2:u2 = " << m1 << ":" << u1 << ":" << m2 << ":" << u2 << endl;
      cout << "gamma(ah1_[i] + bh1_[i]) " << gamma(ah1_[i] + bh1_[i]) << endl;
      cout << "gamma(m1 + u1 + 1.0F) " << gamma(m1 + u1 + 1.0F) << endl;
      cout << "- gamma(ah1_[i] + bh1_[i] + m1 + u1) " << - gamma(ah1_[i] + bh1_[i] + m1 + u1) << endl;
      cout << "gamma(ah1_[i] + m1) " << gamma(ah1_[i] + m1) << endl;
      cout << "- gamma(ah1_[i]) " << - gamma(ah1_[i]) << endl;
      cout << "- gamma(m1 + 1.0F) " << - gamma(m1 + 1.0F) << endl;
      cout << "gamma(bh1_[i] + u1) " << gamma(bh1_[i] + u1) << endl;
      cout << "- gamma(bh1_[i]) " << - gamma(bh1_[i]) << endl;
      cout << "- gamma(u1 + 1.0F) " << - gamma(u1 + 1.0F) << endl;
      */
      ValueType vh = gamma(ah1_[i] + bh1_[i]) + gamma(m1 + u1 + 1.0F) - gamma(ah1_[i] + bh1_[i] + m1 + u1)
	+ gamma(ah1_[i] + m1) - gamma(ah1_[i]) - gamma(m1 + 1.0F) 
	+ gamma(bh1_[i] + u1) - gamma(bh1_[i]) - gamma(u1 + 1.0F);
      ValueType vl = gamma(al1_[i] + bl1_[i]) + gamma(m1 + u1 + 1.0F) - gamma(al1_[i] + bl1_[i] + m1 + u1)
	+ gamma(al1_[i] + m1) - gamma(al1_[i]) - gamma(m1 + 1.0F) 
	+ gamma(bl1_[i] + u1) - gamma(bl1_[i]) - gamma(u1 + 1.0F);
      nmrth[i] = wh1_[i] * Exp(vh);
      nmrtl[i] = wl1_[i] * Exp(vl);
    }
    else {
      nmrth[i] = wh1_[i];
      nmrtl[i] = wl1_[i];
    }
    dnmt += nmrth[i];
    dnmt += nmrtl[i];
  }
  for (uint i=0; i!=wh1_.size(); ++i) {
    pstrh1[i] = nmrth[i] / dnmt;
    pstrl1[i] = nmrtl[i] / dnmt;
  }

  nmrth.clear();
  nmrth.resize(wh2_.size());
  nmrtl.clear();
  nmrtl.resize(wh2_.size());
  dnmt = 0.0;
  for (uint i=0; i!=wh2_.size(); ++i) {
    if (posterior) {
      ValueType vh = gamma(ah2_[i] + bh2_[i]) + gamma(m2 + u2 + 1.0F) - gamma(ah2_[i] + bh2_[i] + m2 + u2)
	+ gamma(ah2_[i] + m2) - gamma(ah2_[i]) - gamma(m2 + 1.0F) 
	+ gamma(bh2_[i] + u2) - gamma(bh2_[i]) - gamma(u2 + 1.0F);
      ValueType vl = gamma(al2_[i] + bl2_[i]) + gamma(m2 + u2 + 1.0F) - gamma(al2_[i] + bl2_[i] + m2 + u2)
	+ gamma(al2_[i] + m2) - gamma(al2_[i]) - gamma(m2 + 1.0F) 
	+ gamma(bl2_[i] + u2) - gamma(bl2_[i]) - gamma(u2 + 1.0F);
      nmrth[i] = wh2_[i] * Exp(vh);
      nmrtl[i] = wl2_[i] * Exp(vl);
    }
    else {
      nmrth[i] = wh2_[i];
      nmrtl[i] = wl2_[i];
    }
    dnmt += nmrth[i];
    dnmt += nmrtl[i];  
  }
  for (uint i=0; i!=wh2_.size(); ++i) {
    pstrh2[i] = nmrth[i] / dnmt;
    pstrl2[i] = nmrtl[i] / dnmt;
  }

  ValueType zval = NEG_INF;
  ValueType val = NEG_INF;
  for (uint i=0; i!=wh1_.size(); ++i) {
    for (uint j=0; j!=wh2_.size(); ++j) {
      ValueType zvhh, zvhl, zvlh, zvll;
      ValueType vhh, vhl, vlh, vll;
      if (posterior) {
	zvhh = Log(pstrh1[i]) + Log(pstrh2[j]) 
	  + beta(m1+m2+ah1_[i]+ah2_[j]-1.0F, u1+u2+bh1_[i]+bh2_[j]-1.0F)
	  - beta(m1+ah1_[i], u1+bh1_[i]) - beta(m2+ah2_[j], u2+bh2_[j]);
	zvhl = Log(pstrh1[i]) + Log(pstrl2[j]) 
	  + beta(m1+m2+ah1_[i]+al2_[j]-1.0F, u1+u2+bh1_[i]+bl2_[j]-1.0F)
	  - beta(m1+ah1_[i], u1+bh1_[i]) - beta(m2+al2_[j], u2+bl2_[j]);
	zvlh = Log(pstrl1[i]) + Log(pstrh2[j]) 
	  + beta(m1+m2+al1_[i]+ah2_[j]-1.0F, u1+u2+bl1_[i]+bh2_[j]-1.0F)
	  - beta(m1+al1_[i], u1+bl1_[i]) - beta(m2+ah2_[j], u2+bh2_[j]);
	zvll = Log(pstrl1[i]) + Log(pstrl2[j]) 
	  + beta(m1+m2+al1_[i]+al2_[j]-1.0F, u1+u2+bl1_[i]+bl2_[j]-1.0F)
	  - beta(m1+al1_[i], u1+bl1_[i]) - beta(m2+al2_[j], u2+bl2_[j]);
	vhh = Log(pstrh1[i]) + Log(pstrh2[j]) 
	  + beta(m1+m2+m1+m2+ah1_[i]+ah2_[j]-1.0F, u1+u2+u1+u2+bh1_[i]+bh2_[j]-1.0F)
	  - beta(m1+ah1_[i], u1+bh1_[i]) - beta(m2+ah2_[j], u2+bh2_[j]);
	vhl = Log(pstrh1[i]) + Log(pstrl2[j]) 
	  + beta(m1+m2+m1+m2+ah1_[i]+al2_[j]-1.0F, u1+u2+u1+u2+bh1_[i]+bl2_[j]-1.0F)
	  - beta(m1+ah1_[i], u1+bh1_[i]) - beta(m2+al2_[j], u2+bl2_[j]);
	vlh = Log(pstrl1[i]) + Log(pstrh2[j]) 
	  + beta(m1+m2+m1+m2+al1_[i]+ah2_[j]-1.0F, u1+u2+u1+u2+bl1_[i]+bh2_[j]-1.0F)
	  - beta(m1+al1_[i], u1+bl1_[i]) - beta(m2+ah2_[j], u2+bh2_[j]);
	vll = Log(pstrl1[i]) + Log(pstrl2[j]) 
	  + beta(m1+m2+m1+m2+al1_[i]+al2_[j]-1.0F, u1+u2+u1+u2+bl1_[i]+bl2_[j]-1.0F)
	  - beta(m1+al1_[i], u1+bl1_[i]) - beta(m2+al2_[j], u2+bl2_[j]);
	/*
	cout << "Log(pstrl1[i]) " << Log(pstrl1[i]) << endl;
	cout << "Log(pstrl2[j]) " << Log(pstrl2[j]) << endl;
	cout << "beta(m1+m2+m1+m2+al1_[i]+al2_[j]-1.0F, u1+u2+u1+u2+bl1_[i]+bl2_[j]-1.0F) " 
	     << beta(m1+m2+m1+m2+al1_[i]+al2_[j]-1.0F, u1+u2+u1+u2+bl1_[i]+bl2_[j]-1.0F) << endl;
	cout << "- beta(m1+al1_[i], u1+bl1_[i]) " << - beta(m1+al1_[i], u1+bl1_[i]) << endl;
	cout << "- beta(m2+al2_[j], u2+bl2_[j]) " << - beta(m2+al2_[j], u2+bl2_[j]) << endl;
	cout << "vhh vhl vlh vll = " << vhh << " " << vhl << " " << vlh << " " << vll << endl;
	*/
      }
      else {
	zvhh = Log(pstrh1[i]) + Log(pstrh2[j]) 
	  + beta(ah1_[i]+ah2_[j]-1.0F, bh1_[i]+bh2_[j]-1.0F)
	  - beta(ah1_[i], bh1_[i]) - beta(ah2_[j], bh2_[j]);
	zvhl = Log(pstrh1[i]) + Log(pstrl2[j]) 
	  + beta(ah1_[i]+al2_[j]-1.0F, bh1_[i]+bl2_[j]-1.0F)
	  - beta(ah1_[i], bh1_[i]) - beta(al2_[j], bl2_[j]);
	zvlh = Log(pstrl1[i]) + Log(pstrh2[j]) 
	  + beta(al1_[i]+ah2_[j]-1.0F, bl1_[i]+bh2_[j]-1.0F)
	  - beta(al1_[i], bl1_[i]) - beta(ah2_[j], bh2_[j]);
	zvll = Log(pstrl1[i]) + Log(pstrl2[j]) 
	  + beta(al1_[i]+al2_[j]-1.0F, bl1_[i]+bl2_[j]-1.0F)
	  - beta(al1_[i], bl1_[i]) - beta(al2_[j], bl2_[j]);
	vhh = Log(pstrh1[i]) + Log(pstrh2[j]) 
	  + beta(m1+m2+ah1_[i]+ah2_[j]-1.0F, u1+u2+bh1_[i]+bh2_[j]-1.0F)
	  - beta(ah1_[i], bh1_[i]) - beta(ah2_[j], bh2_[j]);
	vhl = Log(pstrh1[i]) + Log(pstrl2[j]) 
	  + beta(m1+m2+ah1_[i]+al2_[j]-1.0F, u1+u2+bh1_[i]+bl2_[j]-1.0F)
	  - beta(ah1_[i], bh1_[i]) - beta(al2_[j], bl2_[j]);
	vlh = Log(pstrl1[i]) + Log(pstrh2[j]) 
	  + beta(m1+m2+al1_[i]+ah2_[j]-1.0F, u1+u2+bl1_[i]+bh2_[j]-1.0F)
	  - beta(al1_[i], bl1_[i]) - beta(ah2_[j], bh2_[j]);
	vll = Log(pstrl1[i]) + Log(pstrl2[j]) 
	  + beta(m1+m2+al1_[i]+al2_[j]-1.0F, u1+u2+bl1_[i]+bl2_[j]-1.0F)
	  - beta(al1_[i], bl1_[i]) - beta(al2_[j], bl2_[j]);
      }
      Fast_LogPlusEquals(zval, zvhh);
      Fast_LogPlusEquals(zval, zvhl);
      Fast_LogPlusEquals(zval, zvlh);
      Fast_LogPlusEquals(zval, zvll);
      Fast_LogPlusEquals(val, vhh);
      Fast_LogPlusEquals(val, vhl);
      Fast_LogPlusEquals(val, vlh);
      Fast_LogPlusEquals(val, vll);
    }
  }

  ValueType p = - zval + comb(m1+u1, m1) + comb(m2+u2, m2) + val;
  return p;
}

GlobalStatistics::ValueType GlobalStatistics::
pup_pseudo_count(float m1, float u1, float m2, float u2) const
{
  const float pseudo = 8.0;
  uint m1i = round(m1);
  uint u1i = round(u1);
  uint m2i = round(m2);
  uint u2i = round(u2);
  ValueType q1up = (ValueType) (m1i + pseudo) / (m1i + pseudo + u1i); 
  ValueType q2up = (ValueType) m2i / (m2i + u2i + pseudo); 
  ValueType pup = log_binom(m1i, u1i, q1up) + log_binom(m2i, u2i, q2up);

  return pup;
}


GlobalStatistics::ValueType GlobalStatistics::
pdown_pseudo_count(float m1, float u1, float m2, float u2) const
{
  const float pseudo = 8.0;
  uint m1i = round(m1);
  uint u1i = round(u1);
  uint m2i = round(m2);
  uint u2i = round(u2);
  ValueType q1down = (ValueType) m1i / (m1i + u1i + pseudo); 
  ValueType q2down = (ValueType) (m2i + pseudo) / (m2i + pseudo + u2i); 
  ValueType pdown = log_binom(m1i, u1i, q1down) + log_binom(m2i, u2i, q2down);

  return pdown;
}

GlobalStatistics::ValueType GlobalStatistics::
pnochange_pseudo_count(float m1, float u1, float m2, float u2) const
{
  uint m1i = round(m1);
  uint u1i = round(u1);
  uint m2i = round(m2);
  uint u2i = round(u2);
  ValueType q0 = (ValueType) (m1i + m2i) / (m1i + u1i + m2i + u2i);
  ValueType pnochange = log_binom(m1i, u1i, q0) + log_binom(m2i, u2i, q0);

  return pnochange;
}

GlobalStatistics::ValueType GlobalStatistics::
gamma(float av) const 
{
  assert(av > 0.0);
  
  uint a = (uint) av;
  if (a == av) return gamma(a);
  else return fast_gamma(av);
}

GlobalStatistics::ValueType GlobalStatistics::
digamma(float av) const 
{
  assert(av > 0.0);

  uint a = (uint) av;
  if (a == av) return digamma(a);
  else return fast_digamma(av);
}

GlobalStatistics::ValueType GlobalStatistics::
beta(float x, float y) const
{
  return gamma(x) + gamma(y) - gamma(x + y);
}

GlobalStatistics::ValueType GlobalStatistics::
comb(float x, float y) const 
{
  return gamma(x + 1.0F) - gamma(y + 1.0F) - gamma(x - y + 1.0F);
}

GlobalStatistics::ValueType GlobalStatistics::
gamma(uint a) const 
{
  if (a < gamma_size_) return gamma_[a];
  else return fast_gamma((double) a);
}

GlobalStatistics::ValueType GlobalStatistics::
digamma(uint a) const 
{
  if (a < digamma_size_) return digamma_[a];
  else return fast_digamma((double) a);
}

GlobalStatistics::ValueType GlobalStatistics::
beta(uint a, uint b) const
{
  return gamma(a) + gamma(b) - gamma(a + b);
}

GlobalStatistics::ValueType GlobalStatistics::
comb(uint a, uint b) const 
{
  return gamma(a + 1) - gamma(b + 1) - gamma(a - b + 1);
}

void GlobalStatistics::
test(bool nobeta)
{
  const uint N = 20;

  vector<vector<ValueType> > EmitUp(N+1, vector<ValueType>(N+1, NEG_INF));
  vector<vector<ValueType> > EmitDown(N+1, vector<ValueType>(N+1, NEG_INF));
  vector<vector<ValueType> > EmitNoCh(N+1, vector<ValueType>(N+1, NEG_INF));
  for (uint m1=0; m1<=N; ++m1) {
    for (uint m2=0; m2<=N; ++m2) {
      EmitUp[m1][m2] = pup(m1, N-m1, m2, N-m2);
      EmitDown[m1][m2] = pdown(m1, N-m1, m2, N-m2);
      EmitNoCh[m1][m2] = pnochange(m1, N-m1, m2, N-m2);
    }
  }

  cout << "EmitUp" << endl;
  for (uint m1=0; m1<=N; ++m1) {
    for (uint m2=0; m2<=N; ++m2) {
      //cout << m1 << ":" << m2 << ":";
      cout << EmitUp[m1][m2];
      if (m2 == N) cout << endl;
      else cout << " ";
    }
  }
  cout << "EmitDown" << endl;
  for (uint m1=0; m1<=N; ++m1) {
    for (uint m2=0; m2<=N; ++m2) {
      //cout << m1 << ":" << m2 << ":";
      cout << EmitDown[m1][m2];
      if (m2 == N) cout << endl;
      else cout << " ";
    }
  }
  cout << "EmitNoCh" << endl;
  for (uint m1=0; m1<=N; ++m1) {
    for (uint m2=0; m2<=N; ++m2) {
      //cout << m1 << ":" << m2 << ":";
      cout << EmitNoCh[m1][m2];
      if (m2 == N) cout << endl;
      else cout << " ";
    }
  }

  cout << "EmitUp / EmitNoCh" << endl;
  for (uint m1=0; m1<=N; ++m1) {
    for (uint m2=0; m2<=N; ++m2) {
      //cout << m1 << ":" << m2 << ":";
      cout << formatval("%.1f", EmitUp[m1][m2] - EmitNoCh[m1][m2]);
      if (m2 == N) cout << endl;
      else cout << " ";
    }
  }
  cout << "EmitUp / (EmitDown + EmitNoCh)" << endl;
  for (uint m1=0; m1<=N; ++m1) {
    for (uint m2=0; m2<=N; ++m2) {
      //cout << m1 << ":" << m2 << ":";
      cout << formatval("%.1f", EmitUp[m1][m2] - Fast_LogAdd(EmitDown[m1][m2], EmitNoCh[m1][m2]));
      if (m2 == N) cout << endl;
      else cout << " ";
    }
  }

  cout << "EmitDown / EmitNoCh" << endl;
  for (uint m1=0; m1<=N; ++m1) {
    for (uint m2=0; m2<=N; ++m2) {
      //cout << m1 << ":" << m2 << ":";
      cout << formatval("%.1f", EmitDown[m1][m2] - EmitNoCh[m1][m2]);
      if (m2 == N) cout << endl;
      else cout << " ";
    }
  }
  cout << "EmitDown / (EmitUp + EmitNoCh)" << endl;
  for (uint m1=0; m1<=N; ++m1) {
    for (uint m2=0; m2<=N; ++m2) {
      //cout << m1 << ":" << m2 << ":";
      cout << formatval("%.1f", EmitDown[m1][m2] - Fast_LogAdd(EmitUp[m1][m2], EmitNoCh[m1][m2]));
      if (m2 == N) cout << endl;
      else cout << " ";
    }
  }

}

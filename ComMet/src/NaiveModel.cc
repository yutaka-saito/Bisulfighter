#include <iostream>
#include <queue>
#include <cassert>

#include "NaiveModel.hh"

using namespace std;

void NaiveModel::
reset_param(const MethylList& met, const GlobalStatistics& gstat)
{
  // InitProb
  InitProb.resize(NState);
  for (uint i=0; i!=NState; ++i)
    InitProb[i] = NEG_INF;
  for (uint i=0; i!=NCPGState; ++i)
    InitProb[i] = Log(1.0 / NCPGState);

  // TransProb
  TransProb.resize(NState, vector<ValueType>(NState, NEG_INF));
  const float Sup =    0.99;
  const float Sdown =  0.99;
  const float Snc =    0.9999;
  const float e =      0.1 * (1.0 - Snc) / 2;
  const float f =      0.1 * (1.0 - Snc) / 2;
  const float g =      0.1;
  const float h =      0.1;
  const float Scom =   Snc;
  const float t =      0.25;
  if (noncpg_) {
    TransProb[CPG_UP][CPG_UP] =              Log(t / 3); 
    TransProb[CPG_UP][CPG_DOWN] =            Log(t / 3); 
    TransProb[CPG_UP][CPG_NOCHANGE] =        Log(t / 3); 
    TransProb[CPG_UP][GAP_UP] =              Log(1.0 - g - t); 
    TransProb[CPG_UP][GAP_NOCHANGE] =        Log(g);
    TransProb[CPG_DOWN][CPG_UP] =            Log(t / 3); 
    TransProb[CPG_DOWN][CPG_DOWN] =          Log(t / 3); 
    TransProb[CPG_DOWN][CPG_NOCHANGE] =      Log(t / 3); 
    TransProb[CPG_DOWN][GAP_DOWN] =          Log(1.0 - h - t); 
    TransProb[CPG_DOWN][GAP_NOCHANGE] =      Log(h);
    TransProb[CPG_NOCHANGE][CPG_UP] =        Log(t / 3);
    TransProb[CPG_NOCHANGE][CPG_DOWN] =      Log(t / 3);
    TransProb[CPG_NOCHANGE][CPG_NOCHANGE] =  Log(t / 3);
    TransProb[CPG_NOCHANGE][GAP_NOCHANGE] =  Log(1.0 - t);
    TransProb[GAP_UP][CPG_UP] =              Log(1.0 - Scom);
    TransProb[GAP_UP][GAP_UP] =              Log(Scom);
    TransProb[GAP_DOWN][CPG_DOWN] =          Log(1.0 - Scom);
    TransProb[GAP_DOWN][GAP_DOWN] =          Log(Scom);
    TransProb[GAP_NOCHANGE][CPG_UP] =        Log(e);
    TransProb[GAP_NOCHANGE][CPG_DOWN] =      Log(f);
    TransProb[GAP_NOCHANGE][CPG_NOCHANGE] =  Log(1.0 - e - f - Scom);
    TransProb[GAP_NOCHANGE][GAP_NOCHANGE] =  Log(Scom);
  }
  else {
    TransProb[CPG_UP][GAP_UP] =              Log(1.0 - g); 
    TransProb[CPG_UP][GAP_NOCHANGE] =        Log(g); 
    TransProb[CPG_DOWN][GAP_DOWN] =          Log(1.0 - h);
    TransProb[CPG_DOWN][GAP_NOCHANGE] =      Log(h);
    TransProb[CPG_NOCHANGE][GAP_NOCHANGE] =  Log(1.0);
    TransProb[GAP_UP][CPG_UP] =              Log(1.0 - Sup);
    TransProb[GAP_UP][GAP_UP] =              Log(Sup);
    TransProb[GAP_DOWN][CPG_DOWN] =          Log(1.0 - Sdown);
    TransProb[GAP_DOWN][GAP_DOWN] =          Log(Sdown);
    TransProb[GAP_NOCHANGE][CPG_UP] =        Log(e);
    TransProb[GAP_NOCHANGE][CPG_DOWN] =      Log(f);
    TransProb[GAP_NOCHANGE][CPG_NOCHANGE] =  Log(1.0 - e - f - Snc);
    TransProb[GAP_NOCHANGE][GAP_NOCHANGE] =  Log(Snc);
  }

  // EmitProb
  EmitProb.resize(NState, vector<ValueType>(dpsize_, NEG_INF));
  for (uint i=0; i!=dpsize_; ++i) 
    for (uint j=NCPGState; j!=NState; ++j) 
      EmitProb[j][i] = Log(1.0);
  for (uint i=0; i!=met.pos_size(); ++i) {
    ValueType pup = gstat.pup(met, i); 
    ValueType pdown = gstat.pdown(met, i); 
    ValueType pnochange = gstat.pnochange(met, i); 

    uint ii = met.pos_[i] - met.pos_[0];
    EmitProb[CPG_UP][ii] = pup;
    EmitProb[CPG_DOWN][ii] = pdown;
    EmitProb[CPG_NOCHANGE][ii] = pnochange;
    EmitProb[GAP_UP][ii] = NEG_INF;
    EmitProb[GAP_DOWN][ii] = NEG_INF;
    EmitProb[GAP_NOCHANGE][ii] = NEG_INF;
  }
  //print_table(EmitProb, true);
}

void NaiveModel::
dbase(std::ofstream& ofs, const MethylList& met) 
{
  assert(flg_ppr_);

  for (uint i=0; i!=met.pos_size(); ++i) {
    double theta1 = 0.0;
    double theta2 = 0.0;
    for (uint r=0; r!=met.mc1_[i].size(); ++r) 
      theta1 += met.mc1_[i][r] / met.nc1_[i][r];
    for (uint r=0; r!=met.mc2_[i].size(); ++r) 
      theta2 += met.mc2_[i][r] / met.nc2_[i][r];
    theta1 /= met.mc1_[i].size();
    theta2 /= met.mc2_[i].size();
    ofs << met.name_ << "\t" << met.pos_[i] 
	<< "\t" << formatval("%.8f", theta1) << "\t" << formatval("%.8f", theta2);
    uint ii = met.pos_[i] - met.pos_[0];
    for (uint j=0; j!=NCPGState; ++j) 
      ofs << "\t" << formatval("%.8e", ppr_[j][ii]);
    ofs << endl; 
  }
}

void NaiveModel::
dregion(std::ofstream& ofs, const MethylList& met, const vector<uint>& path) 
{
  bool up = false, down = false;
  uint start = 0, stop = 0;

  for (uint i=0; i!=dpsize_; ++i) {
    if (path[i]==CPG_UP || path[i]==GAP_UP) {
      if (!up) start = i + met.pos_[0];
      if (down) {
	stop = i + met.pos_[0];
	if (stop - start > 1) 
	  ofs << met.name_ << "\t" << start << "\t" << stop << "\tDOWN\t-1" << endl;
      }
      up = true;
      down = false;
    }
    else if (path[i]==CPG_DOWN || path[i]==GAP_DOWN) {
      if (up) {
	stop = i + met.pos_[0];
	if (stop - start > 1)
	  ofs << met.name_ << "\t" << start << "\t" << stop << "\tUP\t+1" << endl;
      }
      if (!down) start = i + met.pos_[0];
      up = false;
      down = true;
    }
    else if (path[i]==CPG_NOCHANGE || path[i]==GAP_NOCHANGE) {
      if (up) {
	stop = i + met.pos_[0];
	if (stop - start > 1)
	  ofs << met.name_ << "\t" << start << "\t" << stop << "\tUP\t+1" << endl;
      }
      else if (down) {
	stop = i + met.pos_[0];
	if (stop - start > 1)
	  ofs << met.name_ << "\t" << start << "\t" << stop << "\tDOWN\t-1" << endl;
      }
      up = false;
      down = false;
    }
    else {
      assert(false);
    }
  }
} 

void NaiveModel::
dregion_viterbi(std::ofstream& ofs, const MethylList& met) 
{
  assert(flg_vtb_ && flg_trc_vtb_);
  dregion(ofs, met, path_vtb_);
}

void NaiveModel::
dregion_pdecode(std::ofstream& ofs, const MethylList& met) 
{
  assert(flg_pdc_ && flg_trc_pdc_);
  dregion(ofs, met, path_pdc_);
}

// compute log P(DMR|UP) / ( P(DMR|NC) + P(DMR|DOWN) ) 
void NaiveModel::
dregion_ratio(std::ofstream& ofs, const MethylList& met, ValueType thsh) 
{
  vector<ValueType> tbl(met.pos_size()); 
  vector<ValueType> tbl_up(met.pos_size()); 
  vector<ValueType> tbl_down(met.pos_size()); 
  vector<ValueType> tbl_nc(met.pos_size()); 
  //deque<bool> trc_tbl(met.pos_size()); // true if DMR was extended
  //deque<bool> trc_msk(met.pos_size(), false); // true if cells were masked
  vector<bool> trc_tbl(met.pos_size()); // true if DMR was extended
  vector<bool> trc_msk(met.pos_size()); // true if cells were masked

  bool allmsk;

  // UP DMR
  trc_msk.assign(met.pos_size(), false);
  while (true) { // each time find the best path  
    tbl.assign(met.pos_size(), NEG_INF);
    tbl.assign(met.pos_size(), NEG_INF);
    tbl_up.assign(met.pos_size(), NEG_INF);
    tbl_down.assign(met.pos_size(), NEG_INF);
    tbl_nc.assign(met.pos_size(), NEG_INF);
    trc_tbl.assign(met.pos_size(), false); 
    allmsk = true;

    if (! trc_msk[0]) {
      allmsk = false;
      bool cpc = (met.pos_[1] - met.pos_[0] == 1) ? true : false;
      tbl_up[0] = TransProb[GAP_NOCHANGE][CPG_UP] + EmitProb[CPG_UP][0];
      tbl_down[0] = TransProb[GAP_NOCHANGE][CPG_DOWN] + EmitProb[CPG_DOWN][0];
      tbl_nc[0] = TransProb[GAP_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][0];
      if (cpc) tbl[0] = tbl_up[0] + Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN], 
						TransProb[CPG_UP][CPG_NOCHANGE]) 
	- Fast_LogAdd(tbl_down[0] + Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP], 
						TransProb[CPG_DOWN][CPG_NOCHANGE]), 
		      tbl_nc[0] + Fast_LogAdd(TransProb[CPG_NOCHANGE][CPG_UP], 
					      TransProb[CPG_NOCHANGE][CPG_DOWN]));
      else tbl[0] = tbl_up[0] + TransProb[CPG_UP][GAP_NOCHANGE] 
	- Fast_LogAdd(tbl_down[0] + TransProb[CPG_DOWN][GAP_NOCHANGE], 
		      tbl_nc[0] + TransProb[CPG_NOCHANGE][GAP_NOCHANGE]);
    }
    for (uint i=1; i<met.pos_size(); ++i) {
      if (trc_msk[i]) continue;
      else allmsk = false;
      uint ii = met.pos_[i] - met.pos_[0];
      uint d = met.pos_[i] - met.pos_[i-1];
      bool cpc = (i != met.pos_size() - 1 && met.pos_[i+1] - met.pos_[i] == 1) ? true : false;
      ValueType val, val_up, val_down, val_nc; // DMR extended
      ValueType v, v_up, v_down, v_nc; // DMR opened
      if (d == 1) {
	val_up = TransProb[CPG_UP][CPG_UP] + EmitProb[CPG_UP][ii];
	val_down = TransProb[CPG_DOWN][CPG_DOWN] + EmitProb[CPG_DOWN][ii];
	val_nc = TransProb[CPG_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][ii];
	v_up = Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP], 
			   TransProb[CPG_NOCHANGE][CPG_UP]) 
	  + EmitProb[CPG_UP][ii];
	v_down = Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN], 
			     TransProb[CPG_NOCHANGE][CPG_DOWN]) 
	  + EmitProb[CPG_DOWN][ii];
	v_nc = Fast_LogAdd(TransProb[CPG_UP][CPG_NOCHANGE], 
			   TransProb[CPG_DOWN][CPG_NOCHANGE]) 
	  + EmitProb[CPG_NOCHANGE][ii];
      }
      else {
	val_up = TransProb[CPG_UP][GAP_UP] 
	  + (d - 2) *  TransProb[GAP_UP][GAP_UP] 
	  + TransProb[GAP_UP][CPG_UP] + EmitProb[CPG_UP][ii];
	val_down = TransProb[CPG_DOWN][GAP_DOWN] 
	  + (d - 2) *  TransProb[GAP_DOWN][GAP_DOWN] 
	  + TransProb[GAP_DOWN][CPG_DOWN] + EmitProb[CPG_DOWN][ii];
	val_nc = TransProb[CPG_NOCHANGE][GAP_NOCHANGE] 
	  + (d - 2) * TransProb[GAP_NOCHANGE][GAP_NOCHANGE] 
	  + TransProb[GAP_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][ii];
	v_up = TransProb[GAP_NOCHANGE][CPG_UP] + EmitProb[CPG_UP][ii];
	v_down = TransProb[GAP_NOCHANGE][CPG_DOWN] + EmitProb[CPG_DOWN][ii];
	v_nc = TransProb[GAP_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][ii];
      }
      if (cpc) {
	val = tbl_up[i-1] + val_up + Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN], 
						 TransProb[CPG_UP][CPG_NOCHANGE])
	  - Fast_LogAdd(tbl_down[i-1] + val_down + Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP], 
							       TransProb[CPG_DOWN][CPG_NOCHANGE]), 
			tbl_nc[i-1] + val_nc + Fast_LogAdd(TransProb[CPG_NOCHANGE][CPG_UP], 
							   TransProb[CPG_NOCHANGE][CPG_DOWN])); 
	v = v_up + Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN], 
			       TransProb[CPG_UP][CPG_NOCHANGE]) 
	  - Fast_LogAdd(v_down + Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP], 
					     TransProb[CPG_DOWN][CPG_NOCHANGE]), 
			v_nc + Fast_LogAdd(TransProb[CPG_NOCHANGE][CPG_UP], 
					   TransProb[CPG_NOCHANGE][CPG_DOWN]));
      } 
      else {
	val = tbl_up[i-1] + val_up + TransProb[CPG_UP][GAP_NOCHANGE] 
	  - Fast_LogAdd(tbl_down[i-1] + val_down + TransProb[CPG_DOWN][GAP_NOCHANGE], 
			tbl_nc[i-1] + val_nc + TransProb[CPG_NOCHANGE][GAP_NOCHANGE]);
	v = v_up + TransProb[CPG_UP][GAP_NOCHANGE] 
	  - Fast_LogAdd(v_down + TransProb[CPG_DOWN][GAP_NOCHANGE], 
			v_nc + TransProb[CPG_NOCHANGE][GAP_NOCHANGE]);
      }
      if (val > v && (! trc_msk[i-1])) { // extension from masked region is not allowed 
	tbl[i] = val;
	tbl_up[i] = tbl_up[i-1] + val_up;
	tbl_down[i] = tbl_down[i-1] + val_down;
	tbl_nc[i] = tbl_nc[i-1] + val_nc;
	trc_tbl[i] = true;
      }
      else {
	tbl[i] = v;
	tbl_up[i] = v_up;
	tbl_down[i] = v_down;
	tbl_nc[i] = v_nc;
	trc_tbl[i] = false;
      }
    }
    if (allmsk) break;

    uint idx = 0, jdx = 0;
    ValueType max = NEG_INF;
    for (uint i=0; i!=met.pos_size(); ++i) {
      if (tbl[i] > max) {
	max = tbl[i];
	jdx = i;
      }
    }
    if (max < thsh) break;
    for (idx=jdx; trc_tbl[idx]; --idx);
    if (idx != jdx) {
      ofs << met.name_ << "\t" << met.pos_[idx] << "\t" << met.pos_[jdx]+1 
	  << "\tUP\t" << max << "\t" << max / (met.pos_[jdx] - met.pos_[idx] + 1) << endl;
    }
    for (uint i=idx; i<=jdx; ++i) trc_msk[i] = true;
  }

  // DOWN DMR
  trc_msk.assign(met.pos_size(), false);
  while (true) { // each time find the best path  
    tbl.assign(met.pos_size(), NEG_INF);
    tbl_up.assign(met.pos_size(), NEG_INF);
    tbl_down.assign(met.pos_size(), NEG_INF);
    tbl_nc.assign(met.pos_size(), NEG_INF);
    trc_tbl.assign(met.pos_size(), false); 
    allmsk = true;

    if (! trc_msk[0]) {
      allmsk = false;
      bool cpc = (met.pos_[1] - met.pos_[0] == 1) ? true : false;
      tbl_up[0] = TransProb[GAP_NOCHANGE][CPG_UP] + EmitProb[CPG_UP][0];
      tbl_down[0] = TransProb[GAP_NOCHANGE][CPG_DOWN] + EmitProb[CPG_DOWN][0];
      tbl_nc[0] = TransProb[GAP_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][0];
      if (cpc) tbl[0] = tbl_down[0] + Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP], 
						  TransProb[CPG_DOWN][CPG_NOCHANGE]) 
	- Fast_LogAdd(tbl_up[0] + Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN], 
					      TransProb[CPG_UP][CPG_NOCHANGE]), 
		      tbl_nc[0] + Fast_LogAdd(TransProb[CPG_NOCHANGE][CPG_UP], 
					      TransProb[CPG_NOCHANGE][CPG_DOWN]));
      else tbl[0] = tbl_down[0] + TransProb[CPG_DOWN][GAP_NOCHANGE] 
	- Fast_LogAdd(tbl_up[0] + TransProb[CPG_UP][GAP_NOCHANGE], 
		      tbl_nc[0] + TransProb[CPG_NOCHANGE][GAP_NOCHANGE]);
    }
    for (uint i=1; i<met.pos_size(); ++i) {
      if (trc_msk[i]) continue;
      else allmsk = false;
      uint ii = met.pos_[i] - met.pos_[0];
      uint d = met.pos_[i] - met.pos_[i-1];
      bool cpc = (i != met.pos_size() - 1 && met.pos_[i+1] - met.pos_[i] == 1) ? true : false;
      ValueType val, val_up, val_down, val_nc; // DMR extended
      ValueType v, v_up, v_down, v_nc; // DMR opened
      if (d == 1) {
	val_up = TransProb[CPG_UP][CPG_UP] + EmitProb[CPG_UP][ii];
	val_down = TransProb[CPG_DOWN][CPG_DOWN] + EmitProb[CPG_DOWN][ii];
	val_nc = TransProb[CPG_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][ii];
	v_up = Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP], 
			   TransProb[CPG_NOCHANGE][CPG_UP]) 
	  + EmitProb[CPG_UP][ii];
	v_down = Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN], 
			     TransProb[CPG_NOCHANGE][CPG_DOWN]) 
	  + EmitProb[CPG_DOWN][ii];
	v_nc = Fast_LogAdd(TransProb[CPG_UP][CPG_NOCHANGE], 
			   TransProb[CPG_DOWN][CPG_NOCHANGE]) 
	  + EmitProb[CPG_NOCHANGE][ii];
      }
      else {
	val_up = TransProb[CPG_UP][GAP_UP] 
	  + (d - 2) *  TransProb[GAP_UP][GAP_UP] 
	  + TransProb[GAP_UP][CPG_UP] + EmitProb[CPG_UP][ii];
	val_down = TransProb[CPG_DOWN][GAP_DOWN] 
	  + (d - 2) *  TransProb[GAP_DOWN][GAP_DOWN] 
	  + TransProb[GAP_DOWN][CPG_DOWN] + EmitProb[CPG_DOWN][ii];
	val_nc = TransProb[CPG_NOCHANGE][GAP_NOCHANGE] 
	  + (d - 2) * TransProb[GAP_NOCHANGE][GAP_NOCHANGE] 
	  + TransProb[GAP_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][ii];
	v_up = TransProb[GAP_NOCHANGE][CPG_UP] + EmitProb[CPG_UP][ii];
	v_down = TransProb[GAP_NOCHANGE][CPG_DOWN] + EmitProb[CPG_DOWN][ii];
	v_nc = TransProb[GAP_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][ii];	
      }
      if (cpc) {
	val = tbl_down[i-1] + val_down + Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP], 
						     TransProb[CPG_DOWN][CPG_NOCHANGE])
	  - Fast_LogAdd(tbl_up[i-1] + val_up + Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN], 
							   TransProb[CPG_UP][CPG_NOCHANGE]), 
			tbl_nc[i-1] + val_nc + Fast_LogAdd(TransProb[CPG_NOCHANGE][CPG_UP], 
							   TransProb[CPG_NOCHANGE][CPG_DOWN])); 
	v = v_down + Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP], 
				 TransProb[CPG_DOWN][CPG_NOCHANGE]) 
	  - Fast_LogAdd(v_up + Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN], 
					   TransProb[CPG_UP][CPG_NOCHANGE]), 
			v_nc + Fast_LogAdd(TransProb[CPG_NOCHANGE][CPG_UP], 
					   TransProb[CPG_NOCHANGE][CPG_DOWN]));
      }
      else {
	val = tbl_down[i-1] + val_down + TransProb[CPG_DOWN][GAP_NOCHANGE] 
	  - Fast_LogAdd(tbl_up[i-1] + val_up + TransProb[CPG_UP][GAP_NOCHANGE], 
			tbl_nc[i-1] + val_nc + TransProb[CPG_NOCHANGE][GAP_NOCHANGE]);	
	v = v_down + TransProb[CPG_DOWN][GAP_NOCHANGE] 
	  - Fast_LogAdd(v_up + TransProb[CPG_UP][GAP_NOCHANGE], 
			v_nc + TransProb[CPG_NOCHANGE][GAP_NOCHANGE]);
      }      
      if (val > v && (! trc_msk[i-1])) { // extension from masked region is not allowed 
	tbl[i] = val;
	tbl_up[i] = tbl_up[i-1] + val_up;
	tbl_down[i] = tbl_down[i-1] + val_down;
	tbl_nc[i] = tbl_nc[i-1] + val_nc;
	trc_tbl[i] = true;
      }
      else {
	tbl[i] = v;
	tbl_up[i] = v_up;
	tbl_down[i] = v_down;
	tbl_nc[i] = v_nc;
	trc_tbl[i] = false;
      }
    }
    if (allmsk) break;

    uint idx = 0, jdx = 0;
    ValueType max = NEG_INF;
    for (uint i=0; i!=met.pos_size(); ++i) {
      if (tbl[i] > max) {
	max = tbl[i];
	jdx = i;
      }
    }
    if (max < thsh) break;
    for (idx=jdx; trc_tbl[idx]; --idx);
    if (idx != jdx) {
      ofs << met.name_ << "\t" << met.pos_[idx] << "\t" << met.pos_[jdx]+1 
	  << "\tDOWN\t" << max << "\t" << max / (met.pos_[jdx] - met.pos_[idx] + 1) << endl;
    }
    for (uint i=idx; i<=jdx; ++i) trc_msk[i] = true;
  }
}

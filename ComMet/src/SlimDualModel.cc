#include <iostream>
#include <queue>
#include <cassert>

#include "SlimDualModel.hh"

using namespace std;

void SlimDualModel::
reset_param(const MethylList& met, const GlobalStatistics& gstat)
{
  // InitProb
  InitProb.resize(NCPGState, Log(1.0 / NCPGState));

  // TransProb
  TransProb.resize(NState, vector<ValueType>(NState, NEG_INF));
  // intra first module
  const float Sup =    0.99;
  const float Sdown =  0.99;
  const float Snc =    0.9999;
  const float e =      0.1 * (1.0 - Snc) / 2;
  const float f =      0.1 * (1.0 - Snc) / 2;
  const float g =      0.1;
  const float h =      0.1;
  const float Scom =   Snc; 
  const float t =      0.25;
  // intra second module
  const float Sup2 =    0.9;
  const float Sdown2 =  0.9;
  const float Snc2 =    0.99;
  const float Scom2 =   Snc2;
  const float t2 =      0.25;
  // inter first module - second module
  const float u   = 0.1 * (1.0 - Snc) / 2;
  const float v   = 0.1 * (1.0 - Snc) / 2;
  const float ur  = 0.1;
  const float vr  = 0.1;
  const float r   = 0.1;
  if (noncpg_) {
    // from first module
    TransProb[CPG_UP][CPG_UP] =               Log(t / 3); 
    TransProb[CPG_UP][CPG_DOWN] =             Log(t / 3); 
    TransProb[CPG_UP][CPG_NOCHANGE] =         Log(t / 3); 
    TransProb[CPG_UP][GAP_UP] =               Log(1.0 - g - t); 
    TransProb[CPG_UP][GAP_NOCHANGE] =         Log(g);
    TransProb[CPG_DOWN][CPG_UP] =             Log(t / 3); 
    TransProb[CPG_DOWN][CPG_DOWN] =           Log(t / 3); 
    TransProb[CPG_DOWN][CPG_NOCHANGE] =       Log(t / 3); 
    TransProb[CPG_DOWN][GAP_DOWN] =           Log(1.0 - h - t); 
    TransProb[CPG_DOWN][GAP_NOCHANGE] =       Log(h);
    TransProb[CPG_NOCHANGE][CPG_UP] =         Log(t / 3);
    TransProb[CPG_NOCHANGE][CPG_DOWN] =       Log(t / 3);
    TransProb[CPG_NOCHANGE][CPG_NOCHANGE] =   Log(t / 3);
    TransProb[CPG_NOCHANGE][GAP_NOCHANGE] =   Log(1.0 - t);
    TransProb[GAP_UP][CPG_UP] =               Log(1.0 - Scom);
    TransProb[GAP_UP][GAP_UP] =               Log(Sup);
    TransProb[GAP_DOWN][CPG_DOWN] =           Log(1.0 - Scom);
    TransProb[GAP_DOWN][GAP_DOWN] =           Log(Scom);
    TransProb[GAP_NOCHANGE][CPG_UP] =         Log(e);
    TransProb[GAP_NOCHANGE][CPG_DOWN] =       Log(f);
    TransProb[GAP_NOCHANGE][CPG_NOCHANGE] =   Log((1.0 - e - f - Scom - u - v) / 2);
    TransProb[GAP_NOCHANGE][CPG_UP2] =        Log(u);
    TransProb[GAP_NOCHANGE][CPG_DOWN2] =      Log(v);
    TransProb[GAP_NOCHANGE][CPG_NOCHANGE2] =  Log((1.0 - e - f - Scom - u - v) / 2);
    TransProb[GAP_NOCHANGE][GAP_NOCHANGE] =   Log(Scom);
    // from second module
    TransProb[CPG_UP2][CPG_UP2] =              Log(t2 / 3); 
    TransProb[CPG_UP2][CPG_DOWN2] =            Log(t2 / 3); 
    TransProb[CPG_UP2][CPG_NOCHANGE2] =        Log(t2 / 3); 
    TransProb[CPG_UP2][GAP_UP2] =              Log(1.0 - ur - t2); 
    TransProb[CPG_UP2][GAP_NOCHANGE] =         Log(ur); 
    TransProb[CPG_DOWN2][CPG_UP2] =            Log(t2 / 3); 
    TransProb[CPG_DOWN2][CPG_DOWN2] =          Log(t2 / 3); 
    TransProb[CPG_DOWN2][CPG_NOCHANGE2] =      Log(t2 / 3); 
    TransProb[CPG_DOWN2][GAP_DOWN2] =          Log(1.0 - vr - t2);
    TransProb[CPG_DOWN2][GAP_NOCHANGE] =       Log(vr);
    TransProb[CPG_NOCHANGE2][CPG_UP2] =        Log(t2 / 3);
    TransProb[CPG_NOCHANGE2][CPG_DOWN2] =      Log(t2 / 3);
    TransProb[CPG_NOCHANGE2][CPG_NOCHANGE2] =  Log(t2 / 3);
    TransProb[CPG_NOCHANGE2][GAP_NOCHANGE2] =  Log(1.0 - r - t2);
    TransProb[CPG_NOCHANGE2][GAP_NOCHANGE] =   Log(r);
    TransProb[GAP_UP2][CPG_UP2] =              Log(1.0 - Scom2);
    TransProb[GAP_UP2][GAP_UP2] =              Log(Scom2);
    TransProb[GAP_DOWN2][CPG_DOWN2] =          Log(1.0 - Scom2);
    TransProb[GAP_DOWN2][GAP_DOWN2] =          Log(Scom2);
    TransProb[GAP_NOCHANGE2][CPG_NOCHANGE2] =  Log(1.0 - Scom2);
    TransProb[GAP_NOCHANGE2][GAP_NOCHANGE2] =  Log(Scom2);
  }
  else {
    // from first module
    TransProb[CPG_UP][GAP_UP] =               Log(1.0 - g); 
    TransProb[CPG_UP][GAP_NOCHANGE] =         Log(g); 
    TransProb[CPG_DOWN][GAP_DOWN] =           Log(1.0 - h);
    TransProb[CPG_DOWN][GAP_NOCHANGE] =       Log(h);
    TransProb[CPG_NOCHANGE][GAP_NOCHANGE] =   Log(1.0);
    TransProb[GAP_UP][CPG_UP] =               Log(1.0 - Sup);
    TransProb[GAP_UP][GAP_UP] =               Log(Sup);
    TransProb[GAP_DOWN][CPG_DOWN] =           Log(1.0 - Sdown);
    TransProb[GAP_DOWN][GAP_DOWN] =           Log(Sdown);
    TransProb[GAP_NOCHANGE][CPG_UP] =         Log(e);
    TransProb[GAP_NOCHANGE][CPG_DOWN] =       Log(f);
    TransProb[GAP_NOCHANGE][CPG_NOCHANGE] =   Log((1.0 - e - f - Snc - u - v) / 2);
    TransProb[GAP_NOCHANGE][CPG_UP2] =        Log(u);
    TransProb[GAP_NOCHANGE][CPG_DOWN2] =      Log(v);
    TransProb[GAP_NOCHANGE][CPG_NOCHANGE2] =  Log((1.0 - e - f - Snc - u - v) / 2);
    TransProb[GAP_NOCHANGE][GAP_NOCHANGE] =   Log(Snc);
    // from second module
    TransProb[CPG_UP2][GAP_UP2] =              Log(1.0 - ur); 
    TransProb[CPG_UP2][GAP_NOCHANGE] =         Log(ur); 
    TransProb[CPG_DOWN2][GAP_DOWN2] =          Log(1.0 - vr);
    TransProb[CPG_DOWN2][GAP_NOCHANGE] =       Log(vr);
    TransProb[CPG_NOCHANGE2][GAP_NOCHANGE2] =  Log(1.0 -r);
    TransProb[CPG_NOCHANGE2][GAP_NOCHANGE] =   Log(r);
    TransProb[GAP_UP2][CPG_UP2] =              Log(1.0 - Sup2);
    TransProb[GAP_UP2][GAP_UP2] =              Log(Sup2);
    TransProb[GAP_DOWN2][CPG_DOWN2] =          Log(1.0 - Sdown2);
    TransProb[GAP_DOWN2][GAP_DOWN2] =          Log(Sdown2);
    TransProb[GAP_NOCHANGE2][CPG_NOCHANGE2] =  Log(1.0 - Snc2);
    TransProb[GAP_NOCHANGE2][GAP_NOCHANGE2] =  Log(Snc2);
  }

  // EmitProb
  EmitProb.resize(NCPGState, vector<ValueType>(dpsize_, NEG_INF));
  for (uint i=0; i!=dpsize_; ++i) {
    ValueType pup = gstat.pup(met, i); 
    ValueType pdown = gstat.pdown(met, i); 
    ValueType pnochange = gstat.pnochange(met, i); 
 
    EmitProb[CPG_UP][i] = pup;
    EmitProb[CPG_DOWN][i] = pdown;
    EmitProb[CPG_NOCHANGE][i] = pnochange;
    EmitProb[CPG_UP2][i] = pup;
    EmitProb[CPG_DOWN2][i] = pdown;
    EmitProb[CPG_NOCHANGE2][i] = pnochange;
  }
  //print_table(EmitProb, true);

  // TransProbDist
  TransProbDist.resize(NGAPState, vector<ValueType>(dpsize_-1, NEG_INF));
  for (uint i=0; i!=dpsize_-1; ++i) {
    for (uint j=NCPGState; j!=NState; ++j) {
      if (Dist[i] == 1) TransProbDist[j-NCPGState][i] = NEG_INF;
      TransProbDist[j-NCPGState][i] = (Dist[i] - 2) * TransProb[j][j];
    }
  }
}

void SlimDualModel::
dbase(std::ofstream& ofs, const MethylList& met) 
{
  assert(flg_ppr_);

  for (uint i=0; i!=dpsize_; ++i) {
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
    for (uint j=0; j!=NCPGState/2; ++j) 
      ofs << "\t" << formatval("%.8e", ppr_[j][i] + ppr_[j+NCPGState/2][i]);
    ofs << endl; 
  }
}

void SlimDualModel::
dregion(std::ofstream& ofs, const MethylList& met, const vector<uint>& path) 
{
  bool up = false, down = false;
  uint start = 0, stop = 0;

  for (uint i=0; i!=dpsize_; ++i) {
    if (path[i]==CPG_UP || path[i]==GAP_UP || 
	path[i]==CPG_UP2 || path[i]==GAP_UP2) {
      if (!up) start = i + met.pos_[0];
      if (down) {
	stop = i + met.pos_[0];
	if (stop - start > 1) 
	  ofs << met.name_ << "\t" << start << "\t" << stop << "\tDOWN\t-1" << endl;
      }
      up = true;
      down = false;
    }
    else if (path[i]==CPG_DOWN || path[i]==GAP_DOWN || 
	     path[i]==CPG_DOWN2 || path[i]==GAP_DOWN2) {
      if (up) {
	stop = i + met.pos_[0];
	if (stop - start > 1)
	  ofs << met.name_ << "\t" << start << "\t" << stop << "\tUP\t+1" << endl;
      }
      if (!down) start = i + met.pos_[0];
      up = false;
      down = true;
    }
    else if (path[i]==CPG_NOCHANGE || path[i]==GAP_NOCHANGE || 
	     path[i]==CPG_NOCHANGE2 || path[i]==GAP_NOCHANGE2) {
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

void SlimDualModel::
dregion_viterbi(std::ofstream& ofs, const MethylList& met) 
{
  assert(flg_vtb_ && flg_trc_vtb_);
  dregion(ofs, met, path_vtb_);
}

void SlimDualModel::
dregion_pdecode(std::ofstream& ofs, const MethylList& met) 
{
  assert(flg_pdc_ && flg_trc_pdc_);
  dregion(ofs, met, path_pdc_);
}


// compute log P(DMR|UP) / ( P(DMR|NC) + P(DMR|DOWN) ) 
void SlimDualModel::
dregion_ratio(std::ofstream& ofs, const MethylList& met, ValueType thsh) 
{
  vector<ValueType> tbl(dpsize_); 
  vector<ValueType> tbl_up(dpsize_); 
  vector<ValueType> tbl_down(dpsize_); 
  vector<ValueType> tbl_nc(dpsize_); 
  vector<ValueType> tbl_up_cgi(dpsize_); 
  vector<ValueType> tbl_down_cgi(dpsize_); 
  vector<ValueType> tbl_nc_cgi(dpsize_); 
  //deque<bool> trc_tbl(dpsize_); // true if DMR was extended
  //deque<bool> trc_msk(dpsize_, false); // true if cells were masked
  vector<bool> trc_tbl(dpsize_); // true if DMR was extended
  vector<bool> trc_msk(dpsize_); // true if cells were masked
  bool allmsk;

  // UP DMR
  trc_msk.assign(dpsize_, false);
  while (true) { // each time find the best path  
    tbl.assign(dpsize_, NEG_INF);
    tbl_up.assign(dpsize_, NEG_INF);
    tbl_down.assign(dpsize_, NEG_INF);
    tbl_nc.assign(dpsize_, NEG_INF);
    tbl_up_cgi.assign(dpsize_, NEG_INF);
    tbl_down_cgi.assign(dpsize_, NEG_INF);
    tbl_nc_cgi.assign(dpsize_, NEG_INF);
    trc_tbl.assign(dpsize_, false); 
    allmsk = true;

    if (! trc_msk[0]) {
      allmsk = false;
      bool cpc = (Dist[0] == 1) ? true : false;
      tbl_up[0] = TransProb[GAP_NOCHANGE][CPG_UP] + EmitProb[CPG_UP][0];
      tbl_down[0] = TransProb[GAP_NOCHANGE][CPG_DOWN] + EmitProb[CPG_DOWN][0];
      tbl_nc[0] = TransProb[GAP_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][0];
      tbl_up_cgi[0] = TransProb[GAP_NOCHANGE][CPG_UP2] + EmitProb[CPG_UP2][0];
      tbl_down_cgi[0] = TransProb[GAP_NOCHANGE][CPG_DOWN2] + EmitProb[CPG_DOWN2][0];
      tbl_nc_cgi[0] = Fast_LogAdd(TransProb[GAP_NOCHANGE][CPG_NOCHANGE2], 
				  TransProb[GAP_NOCHANGE2][CPG_NOCHANGE2]) 
	+ EmitProb[CPG_NOCHANGE2][0];      
      if (cpc) tbl[0] = Fast_LogAdd(tbl_up[0] + Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN], 
							    TransProb[CPG_UP][CPG_NOCHANGE]), // P(DMR|UP) 
				    tbl_up_cgi[0] + Fast_LogAdd(TransProb[CPG_UP2][CPG_DOWN2], 
								TransProb[CPG_UP2][CPG_NOCHANGE2]))
	- Fast_LogAdd(Fast_LogAdd(tbl_down[0] + Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP],
							    TransProb[CPG_DOWN][CPG_NOCHANGE]), // P(DMR|DOWN) 
				  tbl_down_cgi[0] + Fast_LogAdd(TransProb[CPG_DOWN2][CPG_UP2], 
								TransProb[CPG_DOWN2][CPG_NOCHANGE2])),
		      Fast_LogAdd(tbl_nc[0] + Fast_LogAdd(TransProb[CPG_NOCHANGE][CPG_UP],
							  TransProb[CPG_NOCHANGE][CPG_DOWN]), // P(DMR|NC)
				  tbl_nc_cgi[0] + Fast_LogAdd(TransProb[CPG_NOCHANGE2][CPG_UP2], 
							      TransProb[CPG_NOCHANGE2][CPG_DOWN2])));
      else tbl[0] = Fast_LogAdd(tbl_up[0] + TransProb[CPG_UP][GAP_NOCHANGE], // P(DMR|UP) 
				tbl_up_cgi[0] + TransProb[CPG_UP2][GAP_NOCHANGE])
	- Fast_LogAdd(Fast_LogAdd(tbl_down[0] + TransProb[CPG_DOWN][GAP_NOCHANGE], // P(DMR|DOWN) 
				  tbl_down_cgi[0] + TransProb[CPG_DOWN2][GAP_NOCHANGE]),
		      Fast_LogAdd(tbl_nc[0] + TransProb[CPG_NOCHANGE][GAP_NOCHANGE], // P(DMR|NC)
				  tbl_nc_cgi[0] + Fast_LogAdd(TransProb[CPG_NOCHANGE2][GAP_NOCHANGE], 
							      TransProb[CPG_NOCHANGE2][GAP_NOCHANGE2])));
    }
    for (uint i=1; i<dpsize_; ++i) {
      if (trc_msk[i]) continue;
      else allmsk = false;
      uint d = Dist[i-1];
      bool cpc = (i != dpsize_ - 1 && Dist[i] == 1) ? true : false;
      ValueType val, val_up, val_down, val_nc, val_up_cgi, val_down_cgi, val_nc_cgi; // DMR extended
      ValueType v, v_up, v_down, v_nc, v_up_cgi, v_down_cgi, v_nc_cgi; // DMR opened
      if (d == 1) {
	val_up = TransProb[CPG_UP][CPG_UP] + EmitProb[CPG_UP][i] ;
	val_down = TransProb[CPG_DOWN][CPG_DOWN] + EmitProb[CPG_DOWN][i];
	val_nc = TransProb[CPG_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][i];
	val_up_cgi = TransProb[CPG_UP2][CPG_UP2] + EmitProb[CPG_UP2][i] ;
	val_down_cgi = TransProb[CPG_DOWN2][CPG_DOWN2] + EmitProb[CPG_DOWN2][i];
	val_nc_cgi = TransProb[CPG_NOCHANGE2][CPG_NOCHANGE2] + EmitProb[CPG_NOCHANGE2][i];
	v_up = Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP], 
			   TransProb[CPG_NOCHANGE][CPG_UP]) 
	  + EmitProb[CPG_UP][i];
	v_down = Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN], 
			     TransProb[CPG_NOCHANGE][CPG_DOWN]) 
	  + EmitProb[CPG_DOWN][i];
	v_nc = Fast_LogAdd(TransProb[CPG_UP][CPG_NOCHANGE], 
			   TransProb[CPG_DOWN][CPG_NOCHANGE]) 
	  + EmitProb[CPG_NOCHANGE][i];
	v_up_cgi = Fast_LogAdd(TransProb[CPG_DOWN2][CPG_UP2], 
			       TransProb[CPG_NOCHANGE2][CPG_UP2]) 
	  + EmitProb[CPG_UP2][i];
	v_down_cgi = Fast_LogAdd(TransProb[CPG_UP2][CPG_DOWN2], 
				 TransProb[CPG_NOCHANGE2][CPG_DOWN2]) 
	  + EmitProb[CPG_DOWN2][i];
	v_nc_cgi = Fast_LogAdd(TransProb[CPG_UP2][CPG_NOCHANGE2], 
			       TransProb[CPG_DOWN2][CPG_NOCHANGE2]) 
	  + EmitProb[CPG_NOCHANGE2][i];
      }
      else {
	val_up = TransProb[CPG_UP][GAP_UP] 
	  + (d - 2) *  TransProb[GAP_UP][GAP_UP] 
	  + TransProb[GAP_UP][CPG_UP] + EmitProb[CPG_UP][i] ;
	val_down = TransProb[CPG_DOWN][GAP_DOWN] 
	  + (d - 2) *  TransProb[GAP_DOWN][GAP_DOWN] 
	  + TransProb[GAP_DOWN][CPG_DOWN] + EmitProb[CPG_DOWN][i];
	val_nc = TransProb[CPG_NOCHANGE][GAP_NOCHANGE] 
	  + (d - 2) * TransProb[GAP_NOCHANGE][GAP_NOCHANGE] 
	  + TransProb[GAP_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][i];
	val_up_cgi = TransProb[CPG_UP2][GAP_UP2] 
	  + (d - 2) *  TransProb[GAP_UP2][GAP_UP2] 
	  + TransProb[GAP_UP2][CPG_UP2] + EmitProb[CPG_UP2][i] ;
	val_down_cgi = TransProb[CPG_DOWN2][GAP_DOWN2] 
	  + (d - 2) *  TransProb[GAP_DOWN2][GAP_DOWN2] 
	  + TransProb[GAP_DOWN2][CPG_DOWN2] + EmitProb[CPG_DOWN2][i];
	val_nc_cgi = TransProb[CPG_NOCHANGE2][GAP_NOCHANGE2] 
	  + (d - 2) * TransProb[GAP_NOCHANGE2][GAP_NOCHANGE2] 
	  + TransProb[GAP_NOCHANGE2][CPG_NOCHANGE2] + EmitProb[CPG_NOCHANGE2][i];
	v_up = TransProb[GAP_NOCHANGE][CPG_UP] + EmitProb[CPG_UP][i];
	v_down = TransProb[GAP_NOCHANGE][CPG_DOWN] + EmitProb[CPG_DOWN][i];
	v_nc = TransProb[GAP_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][i];
	v_up_cgi = TransProb[GAP_NOCHANGE][CPG_UP2] + EmitProb[CPG_UP2][i];
	v_down_cgi = TransProb[GAP_NOCHANGE][CPG_DOWN2] + EmitProb[CPG_DOWN2][i];
	v_nc_cgi = Fast_LogAdd(TransProb[GAP_NOCHANGE][CPG_NOCHANGE2], 
			       TransProb[GAP_NOCHANGE2][CPG_NOCHANGE2]) 
	  + EmitProb[CPG_NOCHANGE2][i];
      }
      if (cpc) {
	val = Fast_LogAdd(tbl_up[i-1] + val_up + Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN],
							     TransProb[CPG_UP][CPG_NOCHANGE]), // P(DMR|UP) 
			  tbl_up_cgi[i-1] + val_up_cgi + Fast_LogAdd(TransProb[CPG_UP2][CPG_DOWN2], 
								     TransProb[CPG_UP2][CPG_NOCHANGE2]))
	  - Fast_LogAdd(Fast_LogAdd(tbl_down[i-1] + val_down + Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP],
									   TransProb[CPG_DOWN][CPG_NOCHANGE]), // P(DMR|DOWN) 
				    tbl_down_cgi[i-1] + val_down_cgi + Fast_LogAdd(TransProb[CPG_DOWN2][CPG_UP2], 
										   TransProb[CPG_DOWN2][CPG_NOCHANGE2])),
			Fast_LogAdd(tbl_nc[i-1] + val_nc + Fast_LogAdd(TransProb[CPG_NOCHANGE][CPG_UP],
								       TransProb[CPG_NOCHANGE][CPG_DOWN]), // P(DMR|NC)
				    tbl_nc_cgi[i-1] + val_nc_cgi + Fast_LogAdd(TransProb[CPG_NOCHANGE2][CPG_UP2], 
									       TransProb[CPG_NOCHANGE2][CPG_DOWN2])));	
	v = Fast_LogAdd(v_up + Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN],
					   TransProb[CPG_UP][CPG_NOCHANGE]), // P(DMR|UP) 
			v_up_cgi + Fast_LogAdd(TransProb[CPG_UP2][CPG_DOWN2], 
					       TransProb[CPG_UP2][CPG_NOCHANGE2]))
	  - Fast_LogAdd(Fast_LogAdd(v_down + Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP],
							 TransProb[CPG_DOWN][CPG_NOCHANGE]), // P(DMR|DOWN) 
				    v_down_cgi + Fast_LogAdd(TransProb[CPG_DOWN2][CPG_UP2], 
							     TransProb[CPG_DOWN2][CPG_NOCHANGE2])),
			Fast_LogAdd(v_nc + Fast_LogAdd(TransProb[CPG_NOCHANGE][CPG_UP],
						       TransProb[CPG_NOCHANGE][CPG_DOWN]), // P(DMR|NC)
				    v_nc_cgi + Fast_LogAdd(TransProb[CPG_NOCHANGE2][CPG_UP2], 
							   TransProb[CPG_NOCHANGE2][CPG_DOWN2])));
      }
      else {
	val = Fast_LogAdd(tbl_up[i-1] + val_up + TransProb[CPG_UP][GAP_NOCHANGE], // P(DMR|UP) 
			  tbl_up_cgi[i-1] + val_up_cgi + TransProb[CPG_UP2][GAP_NOCHANGE])
	  - Fast_LogAdd(Fast_LogAdd(tbl_down[i-1] + val_down + TransProb[CPG_DOWN][GAP_NOCHANGE], // P(DMR|DOWN) 
				    tbl_down_cgi[i-1] + val_down_cgi + TransProb[CPG_DOWN2][GAP_NOCHANGE]),
			Fast_LogAdd(tbl_nc[i-1] + val_nc + TransProb[CPG_NOCHANGE][GAP_NOCHANGE], // P(DMR|NC)
				    tbl_nc_cgi[i-1] + val_nc_cgi + Fast_LogAdd(TransProb[CPG_NOCHANGE2][GAP_NOCHANGE], 
									       TransProb[CPG_NOCHANGE2][GAP_NOCHANGE2])));	
	v = Fast_LogAdd(v_up + TransProb[CPG_UP][GAP_NOCHANGE], // P(DMR|UP) 
			v_up_cgi + TransProb[CPG_UP2][GAP_NOCHANGE])
	  - Fast_LogAdd(Fast_LogAdd(v_down + TransProb[CPG_DOWN][GAP_NOCHANGE], // P(DMR|DOWN) 
				    v_down_cgi + TransProb[CPG_DOWN2][GAP_NOCHANGE]),
			Fast_LogAdd(v_nc + TransProb[CPG_NOCHANGE][GAP_NOCHANGE], // P(DMR|NC)
				    v_nc_cgi + Fast_LogAdd(TransProb[CPG_NOCHANGE2][GAP_NOCHANGE], 
							   TransProb[CPG_NOCHANGE2][GAP_NOCHANGE2])));	
      }
      if (val > v && (! trc_msk[i-1])) { // extension from masked region is not allowed 
	tbl[i] = val;
	tbl_up[i] = tbl_up[i-1] + val_up;
	tbl_down[i] = tbl_down[i-1] + val_down;
	tbl_nc[i] = tbl_nc[i-1] + val_nc;
	tbl_up_cgi[i] = tbl_up_cgi[i-1] + val_up_cgi;
	tbl_down_cgi[i] = tbl_down_cgi[i-1] + val_down_cgi;
	tbl_nc_cgi[i] = tbl_nc_cgi[i-1] + val_nc_cgi;
	trc_tbl[i] = true;
      }
      else {
	tbl[i] = v;
	tbl_up[i] = v_up;
	tbl_down[i] = v_down;
	tbl_nc[i] = v_nc;
	tbl_up_cgi[i] = v_up_cgi;
	tbl_down_cgi[i] = v_down_cgi;
	tbl_nc_cgi[i] = v_nc_cgi;
	trc_tbl[i] = false;
      }
    }
    if (allmsk) break;

    uint idx = 0, jdx = 0;
    ValueType max = NEG_INF;
    for (uint i=0; i!=dpsize_; ++i) {
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
  trc_msk.assign(dpsize_, false);
  while (true) { // each time find the best path  
    tbl.assign(dpsize_, NEG_INF);
    tbl_up.assign(dpsize_, NEG_INF);
    tbl_down.assign(dpsize_, NEG_INF);
    tbl_nc.assign(dpsize_, NEG_INF);
    tbl_up_cgi.assign(dpsize_, NEG_INF);
    tbl_down_cgi.assign(dpsize_, NEG_INF);
    tbl_nc_cgi.assign(dpsize_, NEG_INF);
    trc_tbl.assign(dpsize_, false); 
    allmsk = true;

    if (! trc_msk[0]) {
      allmsk = false;
      bool cpc = (Dist[0] == 1) ? true : false;
      tbl_up[0] = TransProb[GAP_NOCHANGE][CPG_UP] + EmitProb[CPG_UP][0];
      tbl_down[0] = TransProb[GAP_NOCHANGE][CPG_DOWN] + EmitProb[CPG_DOWN][0];
      tbl_nc[0] = TransProb[GAP_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][0];
      tbl_up_cgi[0] = TransProb[GAP_NOCHANGE][CPG_UP2] + EmitProb[CPG_UP2][0];
      tbl_down_cgi[0] = TransProb[GAP_NOCHANGE][CPG_DOWN2] + EmitProb[CPG_DOWN2][0];
      tbl_nc_cgi[0] = Fast_LogAdd(TransProb[GAP_NOCHANGE][CPG_NOCHANGE2], 
				  TransProb[GAP_NOCHANGE2][CPG_NOCHANGE2]) 
	+ EmitProb[CPG_NOCHANGE2][0];
      if (cpc) tbl[0] = Fast_LogAdd(tbl_down[0] + Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP], 
							      TransProb[CPG_DOWN][CPG_NOCHANGE]), // P(DMR|DOWN) 
				    tbl_down_cgi[0] + Fast_LogAdd(TransProb[CPG_DOWN2][CPG_UP2], 
								  TransProb[CPG_DOWN2][CPG_NOCHANGE2]))
	- Fast_LogAdd(Fast_LogAdd(tbl_up[0] + Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN],
							  TransProb[CPG_UP][CPG_NOCHANGE]), // P(DMR|UP) 
				  tbl_up_cgi[0] + Fast_LogAdd(TransProb[CPG_UP2][CPG_DOWN2], 
							      TransProb[CPG_UP2][CPG_NOCHANGE2])),
		      Fast_LogAdd(tbl_nc[0] + Fast_LogAdd(TransProb[CPG_NOCHANGE][CPG_UP],
							  TransProb[CPG_NOCHANGE][CPG_DOWN]), // P(DMR|NC)
				  tbl_nc_cgi[0] + Fast_LogAdd(TransProb[CPG_NOCHANGE2][CPG_UP2], 
							      TransProb[CPG_NOCHANGE2][CPG_DOWN2])));
      else tbl[0] = Fast_LogAdd(tbl_down[0] + TransProb[CPG_DOWN][GAP_NOCHANGE], // P(DMR|DOWN) 
				tbl_down_cgi[0] + TransProb[CPG_DOWN2][GAP_NOCHANGE])
	- Fast_LogAdd(Fast_LogAdd(tbl_up[0] + TransProb[CPG_UP][GAP_NOCHANGE], // P(DMR|UP) 
				  tbl_up_cgi[0] + TransProb[CPG_UP2][GAP_NOCHANGE]),
		      Fast_LogAdd(tbl_nc[0] + TransProb[CPG_NOCHANGE][GAP_NOCHANGE], // P(DMR|NC)
				  tbl_nc_cgi[0] + Fast_LogAdd(TransProb[CPG_NOCHANGE2][GAP_NOCHANGE], 
							      TransProb[CPG_NOCHANGE2][GAP_NOCHANGE2])));
    }
    for (uint i=1; i<dpsize_; ++i) {
      if (trc_msk[i]) continue;
      else allmsk = false;
      uint d = Dist[i-1];
      bool cpc = (i != dpsize_ - 1 && Dist[i] == 1) ? true : false;
      ValueType val, val_up, val_down, val_nc, val_up_cgi, val_down_cgi, val_nc_cgi; // DMR extended
      ValueType v, v_up, v_down, v_nc, v_up_cgi, v_down_cgi, v_nc_cgi; // DMR opened
      if (d == 1) {
	val_up = TransProb[CPG_UP][CPG_UP] + EmitProb[CPG_UP][i] ;
	val_down = TransProb[CPG_DOWN][CPG_DOWN] + EmitProb[CPG_DOWN][i];
	val_nc = TransProb[CPG_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][i];
	val_up_cgi = TransProb[CPG_UP2][CPG_UP2] + EmitProb[CPG_UP2][i] ;
	val_down_cgi = TransProb[CPG_DOWN2][CPG_DOWN2] + EmitProb[CPG_DOWN2][i];
	val_nc_cgi = TransProb[CPG_NOCHANGE2][CPG_NOCHANGE2] + EmitProb[CPG_NOCHANGE2][i];
	v_up = Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP], 
			   TransProb[CPG_NOCHANGE][CPG_UP]) 
	  + EmitProb[CPG_UP][i];
	v_down = Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN], 
			     TransProb[CPG_NOCHANGE][CPG_DOWN]) 
	  + EmitProb[CPG_DOWN][i];
	v_nc = Fast_LogAdd(TransProb[CPG_UP][CPG_NOCHANGE], 
			   TransProb[CPG_DOWN][CPG_NOCHANGE]) 
	  + EmitProb[CPG_NOCHANGE][i];
	v_up_cgi = Fast_LogAdd(TransProb[CPG_DOWN2][CPG_UP2], 
			       TransProb[CPG_NOCHANGE2][CPG_UP2]) 
	  + EmitProb[CPG_UP2][i];
	v_down_cgi = Fast_LogAdd(TransProb[CPG_UP2][CPG_DOWN2], 
				 TransProb[CPG_NOCHANGE2][CPG_DOWN2]) 
	  + EmitProb[CPG_DOWN2][i];
	v_nc_cgi = Fast_LogAdd(TransProb[CPG_UP2][CPG_NOCHANGE2], 
			       TransProb[CPG_DOWN2][CPG_NOCHANGE2]) 
	  + EmitProb[CPG_NOCHANGE2][i];
      }
      else {
	val_up = TransProb[CPG_UP][GAP_UP] 
	  + (d - 2) *  TransProb[GAP_UP][GAP_UP] 
	  + TransProb[GAP_UP][CPG_UP] + EmitProb[CPG_UP][i];
	val_down = TransProb[CPG_DOWN][GAP_DOWN] 
	  + (d - 2) *  TransProb[GAP_DOWN][GAP_DOWN] 
	  + TransProb[GAP_DOWN][CPG_DOWN] + EmitProb[CPG_DOWN][i];
	val_nc = TransProb[CPG_NOCHANGE][GAP_NOCHANGE] 
	  + (d - 2) * TransProb[GAP_NOCHANGE][GAP_NOCHANGE] 
	  + TransProb[GAP_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][i];
	val_up_cgi = TransProb[CPG_UP2][GAP_UP2] 
	  + (d - 2) *  TransProb[GAP_UP2][GAP_UP2] 
	  + TransProb[GAP_UP2][CPG_UP2] + EmitProb[CPG_UP2][i] ;
	val_down_cgi = TransProb[CPG_DOWN2][GAP_DOWN2] 
	  + (d - 2) *  TransProb[GAP_DOWN2][GAP_DOWN2] 
	  + TransProb[GAP_DOWN2][CPG_DOWN2] + EmitProb[CPG_DOWN2][i];
	val_nc_cgi = TransProb[CPG_NOCHANGE2][GAP_NOCHANGE2] 
	  + (d - 2) * TransProb[GAP_NOCHANGE2][GAP_NOCHANGE2] 
	  + TransProb[GAP_NOCHANGE2][CPG_NOCHANGE2] + EmitProb[CPG_NOCHANGE2][i];
	v_up = TransProb[GAP_NOCHANGE][CPG_UP] + EmitProb[CPG_UP][i];
	v_down = TransProb[GAP_NOCHANGE][CPG_DOWN] + EmitProb[CPG_DOWN][i];
	v_nc = TransProb[GAP_NOCHANGE][CPG_NOCHANGE] + EmitProb[CPG_NOCHANGE][i];
	v_up_cgi = TransProb[GAP_NOCHANGE][CPG_UP2] + EmitProb[CPG_UP2][i];
	v_down_cgi = TransProb[GAP_NOCHANGE][CPG_DOWN2]	+ EmitProb[CPG_DOWN2][i];
	v_nc_cgi = Fast_LogAdd(TransProb[GAP_NOCHANGE][CPG_NOCHANGE2], 
			       TransProb[GAP_NOCHANGE2][CPG_NOCHANGE2]) 
	  + EmitProb[CPG_NOCHANGE2][i];
      }
      if (cpc) {
	val = Fast_LogAdd(tbl_down[i-1] + val_down + Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP],
								 TransProb[CPG_DOWN][CPG_NOCHANGE]), // P(DMR|DOWN) 
			  tbl_down_cgi[i-1] + val_down_cgi + Fast_LogAdd(TransProb[CPG_DOWN2][CPG_UP2], 
									 TransProb[CPG_DOWN2][CPG_NOCHANGE2]))
	  - Fast_LogAdd(Fast_LogAdd(tbl_up[i-1] + val_up + Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN],
								       TransProb[CPG_UP][CPG_NOCHANGE]), // P(DMR|UP) 
				    tbl_up_cgi[i-1] + val_up_cgi + Fast_LogAdd(TransProb[CPG_UP2][CPG_DOWN2], 
									       TransProb[CPG_UP2][CPG_NOCHANGE2])),
			Fast_LogAdd(tbl_nc[i-1] + val_nc + Fast_LogAdd(TransProb[CPG_NOCHANGE][CPG_UP],
								       TransProb[CPG_NOCHANGE][CPG_DOWN]), // P(DMR|NC)
				    tbl_nc_cgi[i-1] + val_nc_cgi + Fast_LogAdd(TransProb[CPG_NOCHANGE2][CPG_UP2], 
									       TransProb[CPG_NOCHANGE2][CPG_DOWN2])));	
	v = Fast_LogAdd(v_down + Fast_LogAdd(TransProb[CPG_DOWN][CPG_UP],
					     TransProb[CPG_DOWN][CPG_NOCHANGE]), // P(DMR|DOWN) 
			v_down_cgi + Fast_LogAdd(TransProb[CPG_DOWN2][CPG_UP2], 
						 TransProb[CPG_DOWN2][CPG_NOCHANGE2]))
	  - Fast_LogAdd(Fast_LogAdd(v_up + Fast_LogAdd(TransProb[CPG_UP][CPG_DOWN],
						       TransProb[CPG_UP][CPG_NOCHANGE]), // P(DMR|UP) 
				    v_up_cgi + Fast_LogAdd(TransProb[CPG_UP2][CPG_DOWN2], 
							   TransProb[CPG_UP2][CPG_NOCHANGE2])),
			Fast_LogAdd(v_nc + Fast_LogAdd(TransProb[CPG_NOCHANGE][CPG_UP],
						       TransProb[CPG_NOCHANGE][CPG_DOWN]), // P(DMR|NC)
				    v_nc_cgi + Fast_LogAdd(TransProb[CPG_NOCHANGE2][CPG_UP2], 
							   TransProb[CPG_NOCHANGE2][CPG_DOWN2])));
      }
      else {
	val = Fast_LogAdd(tbl_down[i-1] + val_down + TransProb[CPG_DOWN][GAP_NOCHANGE], // P(DMR|DOWN) 
			  tbl_down_cgi[i-1] + val_down_cgi + TransProb[CPG_DOWN2][GAP_NOCHANGE])
	  - Fast_LogAdd(Fast_LogAdd(tbl_up[i-1] + val_up + TransProb[CPG_UP][GAP_NOCHANGE], // P(DMR|UP) 
				    tbl_up_cgi[i-1] + val_up_cgi + TransProb[CPG_UP2][GAP_NOCHANGE]),
			Fast_LogAdd(tbl_nc[i-1] + val_nc + TransProb[CPG_NOCHANGE][GAP_NOCHANGE], // P(DMR|NC)
				    tbl_nc_cgi[i-1] + val_nc_cgi + Fast_LogAdd(TransProb[CPG_NOCHANGE2][GAP_NOCHANGE], 
									       TransProb[CPG_NOCHANGE2][GAP_NOCHANGE2])));	
	v = Fast_LogAdd(v_down + TransProb[CPG_DOWN][GAP_NOCHANGE], // P(DMR|DOWN) 
			v_down_cgi + TransProb[CPG_DOWN2][GAP_NOCHANGE])
	  - Fast_LogAdd(Fast_LogAdd(v_up + TransProb[CPG_UP][GAP_NOCHANGE], // P(DMR|UP) 
				    v_up_cgi + TransProb[CPG_UP2][GAP_NOCHANGE]),
			Fast_LogAdd(v_nc + TransProb[CPG_NOCHANGE][GAP_NOCHANGE], // P(DMR|NC)
				    v_nc_cgi + Fast_LogAdd(TransProb[CPG_NOCHANGE2][GAP_NOCHANGE], 
							   TransProb[CPG_NOCHANGE2][GAP_NOCHANGE2])));	
      }
      if (val > v && (! trc_msk[i-1])) { // extension from masked region is not allowed 
	tbl[i] = val;
	tbl_up[i] = tbl_up[i-1] + val_up;
	tbl_down[i] = tbl_down[i-1] + val_down;
	tbl_nc[i] = tbl_nc[i-1] + val_nc;
	tbl_up_cgi[i] = tbl_up_cgi[i-1] + val_up_cgi;
	tbl_down_cgi[i] = tbl_down_cgi[i-1] + val_down_cgi;
	tbl_nc_cgi[i] = tbl_nc_cgi[i-1] + val_nc_cgi;
	trc_tbl[i] = true;
      }
      else {
	tbl[i] = v;
	tbl_up[i] = v_up;
	tbl_down[i] = v_down;
	tbl_nc[i] = v_nc;
	tbl_up_cgi[i] = v_up_cgi;
	tbl_down_cgi[i] = v_down_cgi;
	tbl_nc_cgi[i] = v_nc_cgi;
	trc_tbl[i] = false;
      }
    }
    if (allmsk) break;

    uint idx = 0, jdx = 0;
    ValueType max = NEG_INF;
    for (uint i=0; i!=dpsize_; ++i) {
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


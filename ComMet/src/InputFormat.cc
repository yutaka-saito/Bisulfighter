// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#include <iostream>
#include <sstream>

#include "InputFormat.h"

using namespace std;

bool
load_mb(string& name, vector<uint>& pos, 
	vector<vector<uint> >& mc, vector<vector<uint> >& uc, 
	ifstream& ifs, uint dsep)
{
  static string nm, pre_nm;
  static uint p, pre_p;
  static uint m1, u1, m2, u2;
  static bool init = true;

  if (ifs.eof()) return false;

  if (init) {
    pos.resize(0);
    mc.resize(2, vector<uint>(0));
    uc.resize(2, vector<uint>(0));
  }
  else {
    name = nm;
    pre_nm = nm;
    pos.resize(1, p);
    pre_p = p;
    mc.resize(2, vector<uint>(1));
    mc[0][0] = m1;
    mc[1][0] = m2;
    uc.resize(2, vector<uint>(1));
    uc[0][0] = u1;
    uc[1][0] = u2;
    cout << "new sequence: " << name << endl;
  }

  string s;
  while (getline(ifs, s) != NULL) {
    stringstream ss;
    float m1_buf, u1_buf, m2_buf, u2_buf;
    ss << s;
    ss >> nm >> p >> m1_buf >> u1_buf >> m2_buf >> u2_buf; 
    m1 = round(m1_buf);
    u1 = round(u1_buf);
    m2 = round(m2_buf);
    u2 = round(u2_buf);
    if (m1 + u1 == 0 || m2 + u2 == 0) continue;

    if (init) {
      name = nm;
      init = false;
      cout << "new sequence: " << name << endl;
    }
    else if (nm!=pre_nm || p-pre_p>dsep*1000) { 
      return true;
    }
    pre_nm = nm;
    pre_p = p;
    pos.push_back(p);
    mc[0].push_back(m1);
    mc[1].push_back(m2);
    uc[0].push_back(u1);
    uc[1].push_back(u2);
    //cerr << nm << "\t" << p << "\t" << m1 << "\t" << u1 << "\t" << m2 << "\t" << u2 << endl << flush;
  }

  return true;
}

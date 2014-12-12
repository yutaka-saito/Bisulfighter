#include <iostream>
#include <sstream>
#include <cstdlib>

#include "InputFormat.hh"

using namespace std;

void
split(vector<string>& vec, const string &csv) 
{
  stringstream ss(csv);
  string s;

  vec.clear();
  while(getline(ss, s, ',') != NULL) 
    vec.push_back(s);
}

template < typename StringOrFloat > 
void 
join(string &csv, const vector<StringOrFloat>& vec) 
{
  stringstream ss;

  for (uint i=0; i!=vec.size(); ++i) 
    ss << vec[i] << ",";
  ss >> csv;
  if (! csv.empty()) csv.erase(csv.size() - 1);
}

bool 
parse_counts(vector<float>& m, vector<float>& u, 
	     const string& m_csv, const string& u_csv)
{
  m.clear();
  u.clear();

  vector<string> m_vec;
  vector<string> u_vec;
  split(m_vec, m_csv);
  split(u_vec, u_csv);

  bool ok = (m_vec.size() == u_vec.size()) ? true : false;
  uint size = (m_vec.size() < u_vec.size()) ? m_vec.size() : u_vec.size();
  for (uint i=0; i!=size; ++i) {
    if (m_vec[i].empty() || u_vec[i].empty()) {
      ok = false;
      continue;
    }
    float m_v = atof(m_vec[i].c_str());
    float u_v = atof(u_vec[i].c_str());
    if (m_v < 0.0 || u_v < 0.0) { 
      ok = false;
      continue;
    }
    if (round(m_v) == 0 && round(u_v) == 0) {
      ok = false;
      continue;
    }
    m.push_back(m_v);
    u.push_back(u_v);
  }

  if (m.empty()) {
    cout << "\"" << m_csv << " " << u_csv << "\"" << " is parsed as "
	 << "an empty line" << endl;
    return false;
  }
  if (! ok) {
    string m_tmp;
    string u_tmp;
    join(m_tmp, m);
    join(u_tmp, u);
    cout << "\"" << m_csv << " " << u_csv << "\"" << " is parsed as "
	 << "\"" << m_tmp << " " << u_tmp << "\"" << endl;
  }

  return true;
}

bool
load_mb(string& name, vector<uint>& pos, 
	vector<vector<float> >& mc1, 
	vector<vector<float> >& uc1, 
	vector<vector<float> >& mc2, 
	vector<vector<float> >& uc2, 
	ifstream& ifs, uint dsep)
{
  static string nm_ = "", pre_nm_ = "";
  static uint p_ = 0, pre_p_ = 0;
  static vector<float> m1_, u1_, m2_, u2_;
  static bool init_ = true; // loading the first data

  // prepare for the next file after loading all data
  if (ifs.eof()) {
    nm_ = "";
    pre_nm_ = "";
    p_ = 0;
    pre_p_ = 0;
    m1_.clear();
    u1_.clear();
    m2_.clear();
    u2_.clear();
    init_ = true;

    return false;
  }

  if (init_) {
    name = "";
    pos.clear();
    mc1.clear();
    uc1.clear();
    mc2.clear();
    uc2.clear();
  }
  else {
    name = nm_;
    pre_nm_ = nm_;
    pos.clear();
    pos.resize(1, p_);
    pre_p_ = p_;
    mc1.clear();
    uc1.clear();
    mc2.clear();
    uc2.clear();
    mc1.resize(1, m1_);
    uc1.resize(1, u1_);
    mc2.resize(1, m2_);
    uc2.resize(1, u2_);
    cout << "new sequence: " << name << endl;
  }

  string s;
  while (getline(ifs, s) != NULL) {
    stringstream ss;
    string m1_csv, u1_csv, m2_csv, u2_csv;
    ss << s;
    ss >> nm_ >> p_ >> m1_csv >> u1_csv >> m2_csv >> u2_csv; 
    if (! parse_counts(m1_, u1_, m1_csv, u1_csv)) continue;
    if (! parse_counts(m2_, u2_, m2_csv, u2_csv)) continue;    
    if (init_) {
      name = nm_;
      init_ = false;
      cout << "new sequence: " << name << endl;
    }
    else if (nm_ != pre_nm_ || p_ - pre_p_ > dsep) { 
      return true;
    }
    pre_nm_ = nm_;
    pos.push_back(p_);
    pre_p_ = p_;
    mc1.push_back(m1_);
    uc1.push_back(u1_);
    mc2.push_back(m2_);
    uc2.push_back(u2_);
    /*
    string m1_tmp;
    string u1_tmp;
    string m2_tmp;
    string u2_tmp;
    join(m1_tmp, m1_);
    join(u1_tmp, u1_);
    join(m2_tmp, m2_);
    join(u2_tmp, u2_);
    cerr << nm_ << "\t" << p_ << "\t" << m1_tmp << "\t" << u1_tmp << "\t" << m2_tmp << "\t" << u2_tmp << endl << flush;
    */
  }

  return true;
}

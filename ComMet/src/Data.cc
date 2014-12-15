#include <iostream>
#include <cassert>

#include "Data.hh"
#include "InputFormat.hh"

using namespace std;

MethylList::
MethylList(const string& name, const vector<uint>& pos, 
	   const vector<vector<float> >& mc1, 
	   const vector<vector<float> >& uc1, 
	   const vector<vector<float> >& mc2, 
	   const vector<vector<float> >& uc2)
  : name_(name), pos_(pos), mc1_(mc1), uc1_(uc1), mc2_(mc2), uc2_(uc2)
{
  assert(mc1_.size() == pos_.size() && uc1_.size() == pos_.size());
  assert(mc2_.size() == pos_.size() && uc2_.size() == pos_.size());
  for (uint i=0; i!=pos_.size(); ++i) {
    assert(mc1_[i].size() == uc1_[i].size());
    assert(mc2_[i].size() == uc2_[i].size());
  }

  nc1_.clear();
  nc2_.clear();
  nc1_.resize(pos_.size());
  nc2_.resize(pos_.size());
  for (uint i=0; i!=pos_.size(); ++i) {
    nc1_[i].resize(mc1_[i].size());
    nc2_[i].resize(mc2_[i].size());
    for (uint r=0; r!=nc1_[i].size(); ++r) 
      nc1_[i][r] = mc1_[i][r] + uc1_[i][r];
    for (uint r=0; r!=nc2_[i].size(); ++r) 
      nc2_[i][r] = mc2_[i][r] + uc2_[i][r];
  }

  mc1sum_.clear();
  uc1sum_.clear();
  nc1sum_.clear();
  mc2sum_.clear();
  uc2sum_.clear();
  nc2sum_.clear();
  mc1sum_.resize(pos_.size(), 0.0);
  uc1sum_.resize(pos_.size(), 0.0);
  nc1sum_.resize(pos_.size(), 0.0);
  mc2sum_.resize(pos_.size(), 0.0);
  uc2sum_.resize(pos_.size(), 0.0);
  nc2sum_.resize(pos_.size(), 0.0);
  for (uint i=0; i!=pos_.size(); ++i) {
    for (uint r=0; r!=mc1_[i].size(); ++r) {
      mc1sum_[i] += mc1_[i][r];
      uc1sum_[i] += uc1_[i][r];
      nc1sum_[i] += nc1_[i][r];
    }
    for (uint r=0; r!=mc2_[i].size(); ++r) {
      mc2sum_[i] += mc2_[i][r];
      uc2sum_[i] += uc2_[i][r];
      nc2sum_[i] += nc2_[i][r];
    }
  }
  //print();
}

void MethylList::
print()
{
  for (uint i=0; i!=pos_.size(); ++i) {
    cout << name_ << "\t" << pos_[i] << "\t";
    for (uint r=0; r!=mc1_[i].size(); ++r) {
      cout << mc1_[i][r];
      if (r == mc1_[i].size()-1) cout << "\t";
      else cout << ",";
    }
    cout << "sum:" << mc1sum_[i] << "\t";
    for (uint r=0; r!=uc1_[i].size(); ++r) {
      cout << uc1_[i][r];
      if (r == uc1_[i].size()-1) cout << "\t";
      else cout << ",";
    }
    cout << "sum:" << uc1sum_[i] << "\t";
    for (uint r=0; r!=mc2_[i].size(); ++r) {
      cout << mc2_[i][r];
      if (r == mc2_[i].size()-1) cout << "\t";
      else cout << ",";
    }
    cout << "sum:" << mc2sum_[i] << "\t";
    for (uint r=0; r!=uc2_[i].size(); ++r) {
      cout << uc2_[i][r];
      if (r == uc2_[i].size()-1) cout << "\t";
      else cout << ",";
    }
    cout << "sum:" << uc2sum_[i] << "\t";

    cout << endl << flush;
  }
}

DataLoader::Data* DataLoader::
get(uint dsep)
{
  bool res = false;
  string name;
  vector<uint> pos;
  vector<vector<float> > mc1;
  vector<vector<float> > uc1;
  vector<vector<float> > mc2;
  vector<vector<float> > uc2;

  res = load_mb(name, pos, mc1, uc1, mc2, uc2, ifs_, dsep);
  if (res) return new MethylList(name, pos, mc1, uc1, mc2, uc2);
  else return NULL;
}

DataLoaderFactory::Loader* DataLoaderFactory::
get_loader(ifstream& ifs) 
{
  return new Loader(ifs);
}

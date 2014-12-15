#ifndef __INC_DATA_HH__
#define __INC_DATA_HH__

#include <fstream>
#include <string>
#include <vector>

#include "Utility.hh"

struct MethylList 
{
  std::string name_;
  std::vector<uint> pos_;
  std::vector<std::vector<float> > mc1_;
  std::vector<std::vector<float> > uc1_;
  std::vector<std::vector<float> > nc1_;
  std::vector<std::vector<float> > mc2_;
  std::vector<std::vector<float> > uc2_;
  std::vector<std::vector<float> > nc2_;
  std::vector<float> mc1sum_;
  std::vector<float> uc1sum_;
  std::vector<float> nc1sum_;
  std::vector<float> mc2sum_;
  std::vector<float> uc2sum_;
  std::vector<float> nc2sum_;

  MethylList(const std::string& name, 
	     const std::vector<uint>& pos, 
	     const std::vector<std::vector<float> >& mc1, 
	     const std::vector<std::vector<float> >& uc1, 
	     const std::vector<std::vector<float> >& mc2, 
	     const std::vector<std::vector<float> >& uc2);
  ~MethylList() {}

  void print();

  uint pos_size() const
  {
    return pos_.size();
  }

  uint rep_size(uint c, uint i) const
  {
    assert(c == 1 || c == 2);
    assert(i < pos_.size());

    if (c == 1) return mc1_[i].size();
    else return mc2_[i].size();
  }

  void get_pointer(const std::vector<std::vector<float> >*& mp, 
		   const std::vector<std::vector<float> >*& up, 
		   const std::vector<std::vector<float> >*& np, 
		   const std::vector<float>*& msump, 
		   const std::vector<float>*& usump, 
		   const std::vector<float>*& nsump, 
		   uint c) const
  {
    assert(c == 1 || c == 2);

    if (c == 1) {
      mp = &(mc1_);
      up = &(uc1_);
      np = &(nc1_);
      msump = &(mc1sum_);
      usump = &(uc1sum_);
      nsump = &(nc1sum_);
    }
    else {
      mp = &(mc2_);
      up = &(uc2_);
      np = &(nc2_);
      msump = &(mc2sum_);
      usump = &(uc2sum_);
      nsump = &(nc2sum_);
    }    
  }

  void get_pointer(const std::vector<float>*& mp,
		   const std::vector<float>*& up,
		   const std::vector<float>*& np,
		   const float*& msump, 
		   const float*& usump, 
		   const float*& nsump, 
		   uint c, uint i) const
  {
    assert(c == 1 || c == 2);
    assert(i < pos_.size());

    if (c == 1) {
      mp = &(mc1_[i]);
      up = &(uc1_[i]);
      np = &(nc1_[i]);
      msump = &(mc1sum_[i]);
      usump = &(uc1sum_[i]);
      nsump = &(nc1sum_[i]);
    }
    else {
      mp = &(mc2_[i]);
      up = &(uc2_[i]);
      np = &(nc2_[i]);
      msump = &(mc2sum_[i]);
      usump = &(uc2sum_[i]);
      nsump = &(nc2sum_[i]);
    }    
  }

  void get_pointer(const float*& mp, 
		   const float*& up, 
		   const float*& np, 
		   uint c, uint i, uint r) const
  {
    assert(c == 1 || c == 2);
    assert(i < pos_.size());
    
    if (c == 1) {
      assert(r < mc1_[i].size());
      mp = &(mc1_[i][r]);
      up = &(uc1_[i][r]);
      np = &(nc1_[i][r]);
    }
    else {
      assert(r < mc2_[i].size());
      mp = &(mc2_[i][r]);
      up = &(uc2_[i][r]);
      np = &(nc2_[i][r]);
    }
  }
};

class DataLoader
{
public:
  typedef MethylList Data;

public:
  DataLoader(std::ifstream& ifs) : ifs_(ifs) {}
  ~DataLoader() {}

public:
  Data* get(uint dsep);

private:
  std::ifstream& ifs_;
};

class DataLoaderFactory
{
public:
  typedef DataLoader Loader;
  typedef Loader::Data Data;

public:
  DataLoaderFactory() {}
  ~DataLoaderFactory() {}

public:
  Loader* get_loader(std::ifstream& ifs); 
};

#endif

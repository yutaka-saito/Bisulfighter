// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#ifndef __INC_DATA_H__
#define __INC_DATA_H__

#include <fstream>
#include <vector>
#include <string>

#include "Utility.h"
#include "InputFormat.h"

struct MethylList {
public:
  MethylList(const std::string& name, const std::vector<uint>& pos, 
	     const std::vector<std::vector<uint> >& mc, 
	     const std::vector<std::vector<uint> >& uc) 
    : name_(name), pos_(pos), mc_(mc), uc_(uc) {}

public:
  std::string name_;
  std::vector<uint> pos_;
  std::vector<std::vector<uint> > mc_, uc_;
};

class DataLoader
{
public:
  typedef MethylList Data;

public:
  DataLoader(const std::string& filename);

  Data* get(uint dsep);

private:
  std::string filename_;
  std::ifstream ifs_;
};

class DataLoaderFactory
{
public:
  typedef DataLoader Loader;
  typedef Loader::Data Data;

public:
  DataLoaderFactory() {}

  Loader* get_loader(const std::string filename) const
  {
    return new Loader(filename);
  }
};

#endif

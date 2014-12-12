// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#include <sstream>

#include "Data.h"

using namespace std;
 
DataLoader::
DataLoader(const string& filename)
  : filename_(filename), ifs_(filename.c_str())
{
  if (ifs_.fail()) {
    ostringstream os;
    os << filename_ << ": no such file";
    throw os.str().c_str();
  }
}

DataLoader::Data* DataLoader::
get(uint dsep)
{
  bool ret = false;
  string name;
  vector<uint> pos;
  vector<vector<uint> > mc, uc;

  ret = load_mb(name, pos, mc, uc, ifs_, dsep);
  if (ret) 
    return new MethylList(name, pos, mc ,uc);
  else 
    return NULL;
}

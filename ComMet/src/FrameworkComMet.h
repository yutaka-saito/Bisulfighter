// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#ifndef __INC_FRAMEWORK_COMMET_H__
#define __INC_FRAMEWORK_COMMET_H__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "Utility.h"
#include "ProbabilityModel.h"
#include "Data.h"

struct Options
{
  // input
  std::string imc_file;
  // output
  std::string omc_file;
  std::string dmr_file;
  // others
  static const uint nthr = 1;
  uint nitr;
  uint dsep;
  float thsh;
  float alpha;
  bool dual;
  static const bool noslim = false;
  bool verbose;

  Options() : imc_file(), omc_file(), dmr_file() {}

  void add_options(boost::program_options::options_description& opts);
  void parse_extra_args(const std::vector<std::string>& extra_args);
};

template < class LDF, class MDL >
class App
{
public:
  typedef typename LDF::Data Data;

public:
  App(const LDF& ldf, MDL& mdl, const Options& opts)
    : ldf_(ldf), mdl_(mdl), opts_(opts) {}

  bool execute() 
  {
    if (! (opts_.alpha > 2.0)) {
      std::cout << "alpha must be larger than 2" << std::endl;
      return false;
    }

    return identify();
  }

private:
  bool identify() 
  {
    bool res = false;

    progress("load methylated bases");
    std::vector<Data> data;
    res = load_data(data);
    if (!res) return false;

    std::ofstream ofs_omc(opts_.omc_file.c_str());
    std::ofstream ofs_dmr(opts_.dmr_file.c_str());
    if (ofs_omc.fail()) return false;
    if (ofs_dmr.fail()) return false;
    for (uint i=0; i!=data.size(); ++i) {
      progress("construct a probability model");
      res = mdl_.reset(data[i], opts_.alpha);
      if (!res) continue;

      progress("train parameters");
      mdl_.em(opts_.nitr, opts_.verbose);

      progress("identify differentially methylated bases");
      mdl_.fwdbwd();
      mdl_.posterior();
      mdl_.dbase(ofs_omc, data[i]);

      progress("identify differentially methylated regions");
      mdl_.dregion_ratio(ofs_dmr, data[i], opts_.thsh);
    }

    progress("ComMet - finished");

    return true;
  }

  bool load_data(std::vector<Data>& data)
  {
    typename LDF::Loader* loader=ldf_.get_loader(opts_.imc_file);
    if (loader==NULL) return false;

    std::cout << "loading " << opts_.imc_file << std::endl << std::flush;
    while (true) {
      Data* d = loader->get(opts_.dsep);
      if (d==NULL) break;
      data.push_back(*d);
      delete d;
    }
    delete loader;

    return true;
  }

private:
  const LDF& ldf_;
  MDL& mdl_;
  const Options& opts_;
};

#endif	

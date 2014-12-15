#ifndef __INC_FRAMEWORK_COMMET_HH__
#define __INC_FRAMEWORK_COMMET_HH__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "Utility.hh"
#include "GlobalStatistics.hh"

struct Options
{
  // input
  std::string imc_file;
  // output
  std::string omc_file;
  std::string dmr_file;
  // others
  static const uint nthr = 1;
  uint nmix;
  uint nitr;
  uint nsmp;
  uint dsep;
  float thsh;
  bool noncpg;
  bool nobeta;
  bool nodual;
  bool noslim;
  bool verbose;

  Options() {}
  ~Options() {}

  void add_options(boost::program_options::options_description& opts);
  void parse_extra_args(const std::vector<std::string>& extra_args);
};

template < class LDF, class MDL >
class App
{
public:
  typedef typename LDF::Data Data;

public:
  App(LDF& ldf, MDL& mdl, const Options& opts) : ldf_(ldf), mdl_(mdl), opts_(opts) {}
  ~App() {}

  bool execute() 
  {
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

    progress("calculate global statistics");
    GlobalStatistics gstat;
    gstat.reset(data, opts_.nmix, opts_.nobeta, opts_.nitr, opts_.nsmp, opts_.verbose);
    //gstat.test(opts_.nobeta); return false;

    std::ofstream ofs_omc(opts_.omc_file.c_str());
    std::ofstream ofs_dmr(opts_.dmr_file.c_str());
    if (ofs_omc.fail()) return false;
    if (ofs_dmr.fail()) return false;
    for (uint i=0; i!=data.size(); ++i) {
      progress("construct a probability model");
      res = mdl_.reset(data[i], gstat, opts_.noncpg);
      if (!res) continue;

      progress("estimate parameters");
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
    std::cout << "loading " << opts_.imc_file << std::endl << std::flush;
    std::ifstream ifs_imc(opts_.imc_file.c_str());
    if (ifs_imc.fail()) return false;

    typename LDF::Loader* loader = ldf_.get_loader(ifs_imc);
    while (true) {
      Data* d = loader->get(opts_.dsep);
      if (d == NULL) break;
      data.push_back(*d);
      delete d;
    }
    delete loader;

    return true;
  }

private:
  LDF& ldf_;
  MDL& mdl_;
  const Options& opts_;
};

#endif	

#include <string>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>

#include "Utility.hh"
#include "FrameworkComMet.hh"
#include "NaiveModel.hh"
#include "DualModel.hh"
#include "SlimNaiveModel.hh"
#include "SlimDualModel.hh"
#include "Data.hh"

using namespace std;
using namespace boost;
namespace po = boost::program_options;

int 
main(int argc, char** argv)
{
  progress(whoami());

  // parse command line options
  Options opts;
  po::options_description desc("Options");
  desc.add_options()("help", "show this message");
  opts.add_options(desc);

  po::variables_map vm;
  po::parsed_options parsed =
    po::command_line_parser(argc, argv).
    options(desc).allow_unregistered().run();
  vector<string> extra_args =
    collect_unrecognized(parsed.options, po::include_positional);
  vector<po::option>::iterator new_end =
    remove_if(parsed.options.begin(), parsed.options.end(),
	      bind(&po::option::unregistered, _1) );
  parsed.options.erase(new_end, parsed.options.end());
  po::store(parsed, vm);
  po::notify(vm);
  
  if (vm.count("help") || extra_args.size()<3) {
    cout << endl
	 << "Usage:" << endl
	 << argv[0] << " [options] input output1 output2" << endl
	 << endl
	 << desc << endl;
    return 1;
  }

  opts.parse_extra_args(extra_args);

  if (opts.noncpg) {
    cout << "NOTE: You have specified --noncpg option, which is currently only for testing use." << endl
	 << "See README for some tips about detection of DMRs in non-CpG context." << endl;
    opts.dsep = 1000000;
    opts.nodual = true;
  }

  bool res = false;

  typedef DataLoaderFactory LDF;
  LDF ldf;
  try {
    if (opts.nodual) {
      if (opts.noslim) { 
	typedef NaiveModel MDL;
	MDL mdl;
	App<LDF,MDL> app(ldf,mdl,opts);
	res = app.execute();
      }
      else {
	typedef SlimNaiveModel MDL;
	MDL mdl;
	App<LDF,MDL> app(ldf,mdl,opts);
	res = app.execute();
      }
    }
    else {
      if (opts.noslim) {
	typedef DualModel MDL;
	MDL mdl;
	App<LDF,MDL> app(ldf,mdl,opts);
	res = app.execute();
      }
      else { // default mode
	typedef SlimDualModel MDL;
	MDL mdl;
	App<LDF,MDL> app(ldf,mdl,opts);
	res = app.execute();
      }
    }
  }
  catch (const char* str) {
    cout << str << endl;
  }

  return res ? 0 : 1;
}

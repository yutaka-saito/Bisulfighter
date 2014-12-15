// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#include <string>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>

#include "Utility.h"
#include "FrameworkComMet.h"
#include "NaiveModel.h"
#include "CGIModel.h"
#include "SlimNaiveModel.h"
#include "SlimCGIModel.h"
#include "Data.h"

using namespace std;
using namespace boost;
namespace po = boost::program_options;

int 
main(int argc, char** argv)
{
  progress(whoami());

  Options opts;

  // parse command line options
  po::options_description desc("Options");
  desc.add_options()("help,h", "show this message");
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

  typedef DataLoaderFactory LDF;
  LDF ldf;

  bool res = false;
  try {
    if (opts.dual) {
      if (opts.noslim) {
	typedef CGIModel MDL;
	MDL mdl;
	App<LDF,MDL> app(ldf,mdl,opts);
	res = app.execute();
      }
      else {
	typedef SlimCGIModel MDL;
	MDL mdl;
	App<LDF,MDL> app(ldf,mdl,opts);
	res = app.execute();
      }
    }
    else {
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
  } 
  catch (const char* str) {
    std::cout << str << std::endl;
  }

  return res ? 0 : 1;
}

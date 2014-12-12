#include <cassert>

#include "FrameworkComMet.hh"

using namespace std;
using namespace boost;
namespace po = boost::program_options;

void Options::  
add_options(po::options_description& opts)
{
  opts.add_options()
    //("thread",
    // po::value<uint>(&nthr)->default_value(1),
    // "set the number of threads")
    ("mixture",
     po::value<uint>(&nmix)->default_value(1),
     "set the number of beta mixture components")
    ("iteration",
     po::value<uint>(&nitr)->default_value(500),
     "set the number of training iterations")
    ("sample",
     po::value<uint>(&nsmp)->default_value(10000),
     "set the number of training samples")
    ("separate",
     po::value<uint>(&dsep)->default_value(5000),
     "divide procedures when CpGs are separated by this distance (bp)")
    ("threshold",
     po::value<float>(&thsh)->default_value(0.0),
     "set the threshold for log likelihood ratio scores")
    ("noncpg",
     po::value<bool>(&noncpg)->zero_tokens()->default_value(false),
     "use non-CpG context models (testing; see README)")
    ("nobeta",
     po::value<bool>(&nobeta)->zero_tokens()->default_value(false),
     "DEBUG: use pseudocounts instead of beta mixtures")
    ("nodual",
     po::value<bool>(&nodual)->zero_tokens()->default_value(false),
     "DEBUG: use the naive model instead of the dual model")
    ("noslim",
     po::value<bool>(&noslim)->zero_tokens()->default_value(false),
     "DEBUG: use the straightforward impl instead of the slim impl")
    ("verbose",
     po::value<bool>(&verbose)->zero_tokens()->default_value(false),
     "make verbose reports to stdout")
    ;
}

void Options::
parse_extra_args(const vector<string>& extra_args)
{
  assert(extra_args.size() >= 3);
  imc_file = extra_args[0];
  omc_file = extra_args[1];
  dmr_file = extra_args[2];
}

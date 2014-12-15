// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#include <cassert>

#include "FrameworkComMet.h"

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
    ("iteration",
     po::value<uint>(&nitr)->default_value(100),
     "set the number of training iterations")
    ("separate",
     po::value<uint>(&dsep)->default_value(5),
     "divide procedures when CpGs are separated by this distance (Kb)")
    ("threshold",
     po::value<float>(&thsh)->default_value(0.0),
     "set the threshold for log likelihood ratio scores")
    ("alpha",
     po::value<float>(&alpha)->default_value(8.0),
     "set pseudocounts for the emission function")
    ("dual",
     po::value<bool>(&dual)->zero_tokens()->default_value(false),
     "use the dual model instead of the naive model")
    //("noslim",
    // po::value<bool>(&noslim)->zero_tokens()->default_value(false),
    // "use straightforward implementation instead of slim one")
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


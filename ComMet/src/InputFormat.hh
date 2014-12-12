#ifndef  __INC_INPUT_FORMAT_HH__
#define  __INC_INPUT_FORMAT_HH__

#include <fstream>
#include <string>
#include <vector>

#include "Utility.hh"

bool 
load_mb(std::string& name, std::vector<uint>& pos, 
	std::vector<std::vector<float> >& mc1, 
	std::vector<std::vector<float> >& uc1, 
	std::vector<std::vector<float> >& mc2, 
	std::vector<std::vector<float> >& uc2, 
	std::ifstream& ifs, uint dsep);

#endif 


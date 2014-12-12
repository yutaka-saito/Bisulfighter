// ComMet
// by National Institute of Advanced Industrial Science and Technology (AIST)
// is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// http://creativecommons.org/licenses/by-nc-sa/3.0/


#ifndef  __INC_INPUT_FORMAT_H__
#define  __INC_INPUT_FORMAT_H__

#include <fstream>
#include <string>
#include <vector>

#include "Utility.h"

bool 
load_mb(std::string& name, std::vector<uint>& pos, 
	std::vector<std::vector<uint> >& mc, 
	std::vector<std::vector<uint> >& uc, 
	std::ifstream& ifs, uint dsep);

#endif 


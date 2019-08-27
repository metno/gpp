#include "Simple.h"
#include <fstream>
#include <sstream>
#include "../Util.h"
#include <assert.h>
#include <set>
#include <fstream>

ParameterFileSimple::ParameterFileSimple(Parameters iParameters) : ParameterFile(Options()) {
   Location defaultLocation(Util::MV,Util::MV,Util::MV);
   setParameters(iParameters, 0, defaultLocation);

   recomputeTree();
}

std::vector<int> ParameterFileSimple::getTimes() const {
   std::vector<int> times(1,0);
   return times;
}


bool ParameterFileSimple::isValid(std::string iFilename) {
   // TODO
   return true;
}

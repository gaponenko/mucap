// Andrei Gaponenko, 2013

#include "MuCapUtilities/inc/TabulatedFunction.hh"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

// c++11 regex is not yet supported by our gcc
#include <boost/regex.hpp>

namespace mucap {

  //================================================================
  TabulatedFunction::TabulatedFunction(const std::string& dataFileName) {

    std::ifstream infile(dataFileName.c_str());
    if (!(infile.is_open())) {
      std::ostringstream os;
      os<<"TabulatedFunction(): can't read input file \""<<dataFileName<<"\"\n";
      throw std::runtime_error(os.str());
    }

    std::string line;
    const boost::regex COMMENT("^\\s*#");
    while(getline(infile, line)) {
      if(!boost::regex_search(line, COMMENT)) {
        double e{0}, p{0};
        std::istringstream is(line);
        is.exceptions(std::ios_base::badbit|std::ios_base::failbit);
        is>>e>>p;
        table_.emplace_back(e,p);
      }
    }

    std::sort(table_.begin(), table_.end(),
              [](const Point& a, const Point& b) { return a.x < b.x; }
              );
  }

  //================================================================
  double TabulatedFunction::operator()(double x) const {

    const auto pos = 
      std::lower_bound(table_.begin(), table_.end(), x,
                       [](const Point& a, double x) { return a.x < x; }
                       );

    return  ((pos != table_.end()) && (pos + 1 != table_.end())) ?
      interpolate1(x, *pos, *(pos+1)) : 0.;
  }

  //================================================================
  double TabulatedFunction::interpolate1(double x, const Point& p1, const Point& p2) const {
    return p1.y + (x - p1.x)*(p2.y - p1.y)/(p2.x - p1.x);
  }

  //================================================================
}

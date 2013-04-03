// A generic tabulated function.
//
// Andrei Gaponenko, 2013
//

#ifndef MuCapUtilities_TabulatedFunction_hh
#define MuCapUtilities_TabulatedFunction_hh

#include <string>
#include <vector>

namespace mucap {

  class TabulatedFunction {
  public:

    struct Point {
      double x;
      double y;
      Point(double xx, double yy): x(xx), y(yy) {}
    };

    typedef std::vector<Point> Table;
    const Table& table() const { return table_; }

    explicit TabulatedFunction(const std::string& dataFileName);

    // Interpolated value.  Returns 0 for x out of range.
    double operator()(double x) const;

  private:

    Table table_;

    // linear interpolation
    double interpolate1(double x, const Point& p1, const Point& p2) const;
  };

} // end of namespace mu2e

#endif /* MuCapUtilities_TabulatedFunction_hh */

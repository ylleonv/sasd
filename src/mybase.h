#include <Rcpp.h>
// example from:
// http://stackoverflow.com/questions/24317910/rcpp-module-for-inheritance-class

class Base {
public:
  Base(double x_) : x(x_){}
  double x;
  void addOneToX();
} ;
#ifndef MYBASE_h
#define MYBASE_h

#include <Rcpp.h>

class Base {
public:
  Base(double x_) : x(x_){}
  double x;
  void addOneToX();
} ;

#endif
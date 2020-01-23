#include <Rcpp.h>
#include "mybase.h"

class Derived : public Base {
public:
  Derived(int y_) : Base(0.0), y(y_){}
  int y ;
  void addOneToY();
} ;
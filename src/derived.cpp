#include <Rcpp.h>
#include "derived.h"

void Derived::addOneToY(){
  y = y + 1;
}

RCPP_MODULE(DERIVEDMODULE){
  Rcpp::class_<Base>("Base")
  .field("x", &Base::x)
  .method("addOneToX", &Base::addOneToX)
  ;
  Rcpp::class_<Derived>("Derived")
    .derives<Base>("Base")
    .constructor<int>()
    .field("y", &Derived::y)
    .method("addOneToY", &Derived::addOneToY)
  ;
}
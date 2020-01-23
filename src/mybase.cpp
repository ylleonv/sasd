#include <Rcpp.h>
#include "mybase.h"

//simple method to add one
void Base::addOneToX(){
  x = x + 1;
}

RCPP_MODULE(BASEMODULE){
  using namespace Rcpp ;
  class_<Base>("Base")
    .constructor<double>()
    .field("x", &Base::x)
    .method("addOneToX", &Base::addOneToX)
  ;
}
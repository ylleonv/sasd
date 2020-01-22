#ifndef GUARD_LocationClass_h
#define GUARD_LocationClass_h
#include<Rcpp.h>

class Location{
public:
  Location();
  Location(int, int);
  int x;
  int y;
  void print();
};

RCPP_MODULE(mimod){
  using namespace Rcpp ;
  class_<Location>("Location")
    // expose the default constructor
    .constructor()
    .constructor<int,int>("documentation for constructor")
    .field( "x", &Location::x, "documentation for x")
    .field( "y", &Location::y, "documentation for y")
    .method( "print", &Location::print, "documentation for print")
    ;
}
#endif
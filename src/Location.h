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
#endif
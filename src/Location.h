#ifndef LocationClass_h
#define LocationClass_h
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
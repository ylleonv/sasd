#include "Location.h"
#include<Rcpp.h>

//constructors
Location::Location() : x(0), y(0) { }
Location::Location(int xi, int yi) : x(xi), y(yi) { }

//print function
void Location::print(){
  Rcpp::Rcout << "Lorena = " << x << std::endl;
  Rcpp::Rcout << "y = " << y << std::endl;
}

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


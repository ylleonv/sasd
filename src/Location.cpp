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


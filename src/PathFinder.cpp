#include "PathFinder.h"
#include <math.h>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"


PathFinder::PathFinder() {
  
  //const double STEER_LIMIT     = 50; // * pi() / 180; 
  
  
}

PathFinder::~PathFinder() {}


double PathFinder::myfunc(double deg) {
    return deg * M_PI / 180.0;
}

//###########helpers########################

// A function that returns a value between 0 and 1 for x in the 
// range [0, infinity] and -1 to 1 for x in the range [-infinity, infinity].
// Useful for cost functions.

double PathFinder::logistic(double x) {
  return 2.0 / (1 + exp(-x)) - 1.0;
}

/*
def logistic(x):
    """
    A function that returns a value between 0 and 1 for x in the 
    range [0, infinity] and -1 to 1 for x in the range [-infinity, infinity].

    Useful for cost functions.
    """
      return 2.0 / (1 + exp(-x)) - 1.0
*/  
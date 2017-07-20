#include "PathFinder.h"
#include <math.h>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"


PathFinder::PathFinder() {
  
  //const double STEER_LIMIT     = 50; // * pi() / 180; 
  
  
}

PathFinder::~PathFinder() {}
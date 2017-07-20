#ifndef PATHFINDER_H
#define PATHFINDER_H

#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "spline.h" // rbx

class PathFinder 
{
  private:
    
  
  public:
        
    PathFinder();

    virtual ~PathFinder();
        
    const double STEER_LIMIT     = 5;
    
};

#endif /* PATHFINDER_H */
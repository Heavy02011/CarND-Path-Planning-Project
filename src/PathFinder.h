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
        
    // number of sample points on track
    const int N_SAMPLES = 10;
      
    // set max speed
    const double SPEED_LIMIT = 50.0;
    
    // vehicle radius, model vehicle as circle to simplify collision detection
    const double VEHICLE_RADIUS = 1.5; 
      
    // max jerk and acceleration
    const double MAX_JERK = 10.0; // m/s/s/s
    const double MAX_ACCEL= 10.0; // m/s/s
    
    // max jerk and acceleration in one second
    const double EXPECTED_JERK_IN_ONE_SEC = 2.0; // m/s/s
    const double EXPECTED_ACC_IN_ONE_SEC = 1.0; // m/s
      
    
};

#endif /* PATHFINDER_H */

/*
SIGMA_S = [10.0, 4.0, 2.0] # s, s_dot, s_double_dot
SIGMA_D = [1.0, 1.0, 1.0]
SIGMA_T = 2.0
*/
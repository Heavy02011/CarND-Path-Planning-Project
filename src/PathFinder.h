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

#include <cmath>
//#include "Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

struct test_case {
	
		vector<double> start;
		vector<double> end;
		double T;
};

class PathFinder 
{

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
    
    // sigmas for generating perturbed goal points according to lesson 5.30
    //const vector<double> SIGMA_S = {10.0, 4.0, 2.0}; // s, s_dot, s_double_dot
    //const vector<double> SIGMA_D = {1.0, 1.0, 1.0};
    VectorXd SIGMA_S = VectorXd(3);
    VectorXd SIGMA_D = VectorXd(3);
    const double SIGMA_T = 2.0;
  
    // sample function
    double myfunc(double deg);
    
    // Calculate the Jerk Minimizing Trajectory that connects the initial state to the final state in time T.
    vector<double> JMT(vector< double> start, vector <double> end, double T);
    
    // generate straight path
    void generate_path_straight(double car_x, double car_y, double car_yaw, vector<double> *next_x, vector<double> *next_y); //, vector<double> *prev_x, vector<double> *prev_y);
    
    //===== helper functions =====
    double logistic(double x);
      
  private:
      
      
};

#endif /* PATHFINDER_H */
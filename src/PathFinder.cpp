#include "PathFinder.h"
#include <math.h>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
//#include "Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


PathFinder::PathFinder() {
  
  // sigmas for generating perturbed goal points according to lesson 5.30
  SIGMA_S << 10.0, 4.0, 2.0 ;
  SIGMA_D << 1.0, 1.0, 1.0 ;
      
  
}

PathFinder::~PathFinder() {}

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad2(double x) { return x * pi() / 180; }
double rad2deg2(double x) { return x * 180 / pi(); }


double PathFinder::myfunc(double deg) {
    return deg * M_PI / 180.0;
}

vector<double> PathFinder::JMT(vector< double> start, vector <double> end, double T)
{
    /*
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
        length three array.

    T     - The duration, in seconds, over which this maneuver should occur.

    OUTPUT 
    an array of length 6, each value corresponding to a coefficent in the polynomial 
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

    EXAMPLE

    > JMT( [0, 10, 0], [10, 10, 0], 1)
    [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
    */
    
    // prepare matrix A with coefficents
    int w = 3;
    int h = 3;
    MatrixXd A(w,h);
    
    A(0,0) = T*T*T;
    A(0,1) = T*T*T*T;
    A(0,2) = T*T*T*T*T;
    A(1,0) = 3*T*T;
    A(1,1) = 4*T*T*T;
    A(1,2) = 5*T*T*T*T;
    A(2,0) = 6*T;
    A(2,1) = 12*T*T;
    A(2,2) = 20*T*T*T;
    //cout << A << endl;
    
    // prepare coefficents for vector
    double si         = start[0]; // 0
    double si_dot     = start[1]; // 10
    double si_dot_dot = start[2]; // 0
    double sf         = end[0];   // 10
    double sf_dot     = end[1];   // 10
    double sf_dot_dot = end[2];   // 0
    
    VectorXd b(3);
    b(0) = sf - (si + si_dot*T + 0.5*si_dot_dot*T*T); // 0
    //B[1] = sf_dot - (si_dot*T + si_dot_dot*T) # sf_dot - (si_dot + si_ddot * T)
    b(1) = sf_dot - (si_dot + si_dot_dot*T); // 0
    b(2) = sf_dot_dot - si_dot_dot; // 0
    //cout << b << endl;
    
    // solve the systemof lenar equations using eigen library
    VectorXd x = A.colPivHouseholderQr().solve(b);
    //cout << "The solution is:\n" << x << endl;
    
    // compose solution
    double a0 = si;
    double a1 = si_dot;
    double a2 = si_dot_dot / 2.0;
    double a3 = x[0];
    double a4 = x[1];
    double a5 = x[2];
    
    return {a0,a1,a2,a3,a4,a5};
    
}

void PathFinder::generate_path_straight(double car_x, double car_y, double car_yaw, vector<double> *next_x, vector<double> *next_y) {//, vector<double> *prev_x, vector<double> *prev_y) {
  
  // https://stackoverflow.com/questions/13295011/altering-the-values-of-array-elements-from-within-a-function
  
  //vector<double> next_x_vals;
  //vector<double> next_y_vals;
  double my_x;
  double my_y;
  
  double dist_inc = 0.5;
  for(int i = 0; i < 50; i++)
  {
    //next_x_vals.push_back(car_x+(dist_inc*i)*cos(deg2rad2(car_yaw)));
    //next_y_vals.push_back(car_y+(dist_inc*i)*sin(deg2rad2(car_yaw)));
    my_x = car_x+(dist_inc*i)*cos(deg2rad2(car_yaw));
    my_y = car_y+(dist_inc*i)*sin(deg2rad2(car_yaw));
    (*next_x)[i] = my_x;
    (*next_y)[i] = my_y; 
    
  }
  //cout << "next_x/y_vals.size() = " << next_x_vals.size() << " " << next_y_vals.size() << endl;
  cout << "next_x/y.size() = " << next_x->size() << " " << next_y->size() << endl;
    
}

void PathFinder::generate_path_circle(double car_x, double car_y, double car_yaw, vector<double> *next_x, vector<double> *next_y, vector<double> *prev_x, vector<double> *prev_y) {
  // https://stackoverflow.com/questions/13295011/altering-the-values-of-array-elements-from-within-a-function
  
  vector<double> next_x_vals;
  vector<double> next_y_vals;
  
  double pos_x;
  double pos_y;
  double angle;
  int path_size = prev_x->size();
  cout << "circle: path_size = " << path_size << endl;

  for(int i = 0; i < path_size; i++)
  {
      next_x_vals.push_back((*prev_x)[i]);
      next_y_vals.push_back((*prev_y)[i]);
  }

  if(path_size == 0)
  {
      pos_x = car_x;
      pos_y = car_y;
      angle = deg2rad2(car_yaw);
  }
  else
  {
      pos_x = (*prev_x)[path_size-1];
      pos_y = (*prev_y)[path_size-1];

      double pos_x2 = (*prev_x)[path_size-2];
      double pos_y2 = (*prev_y)[path_size-2];
      angle = atan2(pos_y-pos_y2,pos_x-pos_x2);
  }

  double dist_inc = 0.5;
  for(int i = 0; i < 50-path_size; i++)
  {    
      next_x_vals.push_back(pos_x+(dist_inc)*cos(angle+(i+1)*(pi()/100)));
      next_y_vals.push_back(pos_y+(dist_inc)*sin(angle+(i+1)*(pi()/100)));
      pos_x += (dist_inc)*cos(angle+(i+1)*(pi()/100));
      pos_y += (dist_inc)*sin(angle+(i+1)*(pi()/100));
  }
 
  // store final result
 for(int i=0; i< path_size; i++) {
   (*next_x)[i] = next_x_vals[i];
   (*next_y)[i] = next_y_vals[i];
 } 
  
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
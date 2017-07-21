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
  
  //const double STEER_LIMIT     = 50; // * pi() / 180; 
  
  
}

PathFinder::~PathFinder() {}


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
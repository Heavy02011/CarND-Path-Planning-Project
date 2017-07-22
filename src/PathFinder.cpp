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


PathFinder::PathFinder(double x, double y, double d, double s, double v, double a, double yaw) {
 
  // store data of my car here
  this->d = d;
  this->s = s;
  this->v = v;
  this->a = a;
  this->yaw = yaw;
  state = KL;
  max_acceleration = -1;
  
}

PathFinder::~PathFinder() {}

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad2(double x) { return x * pi() / 180; }
double rad2deg2(double x) { return x * 180 / pi(); }


double PathFinder::myfunc(double deg) {
    return deg * M_PI / 180.0;
}

vector<double> PathFinder::JMT(vector<double> start, vector <double> end, double T)
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

vector<states> PathFinder::successor_states(states input)
{

  vector<states> succ_states;
  switch(state)
  {
    case KL: // 0
        succ_states.push_back(KL);   // 0
        succ_states.push_back(PLCL); // 1     
        succ_states.push_back(PLCR); // 2       
        return succ_states;
    case PLCL: // 1
        succ_states.push_back(PLCL); // 1
        succ_states.push_back(KL);   // 0
        succ_states.push_back(LCL);  // 3      
        return succ_states;
    case PLCR: // 2
        succ_states.push_back(PLCR); // 2 
        succ_states.push_back(KL);   // 0
        succ_states.push_back(LCR);  // 4      
        return succ_states;
    case LCL: // 3
        succ_states.push_back(LCL);  // 3
        succ_states.push_back(KL);   // 0   
        return succ_states;
    case LCR: // 4
        succ_states.push_back(LCR);  // 4
        succ_states.push_back(KL);   // 0
        return succ_states;
    default:
        cerr << "PathFinder::successor_states: state not defined -> " << state << endl;
        break;
  }
}


//###########helpers########################

// A function that returns a value between 0 and 1 for x in the 
// range [0, infinity] and -1 to 1 for x in the range [-infinity, infinity].
// Useful for cost functions.
double PathFinder::logistic(double x) {
  return 2.0 / (1 + exp(-x)) - 1.0;
}

// pass a vector x to cout
void PathFinder::output_vector(vector<double> x) {
  for (int i=0; i<x.size(); i++) {
    cout << i << " " << x[i] << endl;  
  } 
  cout << endl;
}

// evaluate a polynomal 5th order determined by coefficients at value x
double PathFinder::evaluate_polynomal(vector<double> coeffcients, double x) {
  double function = 0.0;
  for (int i = 0; i < coeffcients.size(); i++) {
    function += coeffcients[i] * pow(x, i);
  }
  return function;
}

/*
function<int(double)> PathFinder::f(char c) {
  double g(double x) {
    g = 23.0*x;
    return g*g;
  }
}
*/

/*
function<double(vector<double>)> PathFinder::to_equation(vector<double> coefficents) {
  
  double function(vector<double> coefficents, double x) {
    for (int i = 0; i < coeffcients.size(); i++) {
      function += coeffcients[i] * pow(x, i);
    }
  }
  return function;
}
*/

//calculates the derivative of a polynomial and returns the corresponding coefficients. ref: helpers.py in lesson 5.30
vector<double> PathFinder::differentiate(vector<double> coefficients) {
  vector<double> new_coeff;
  for (int i = 1; i < coefficients.size(); i++) {
    new_coeff.push_back( (float(i)+1.0) * coefficients[i] );
  }
  return new_coeff;
}

// calculate distance of othercar to my car
double PathFinder::distance2car(Vehicle othercar) {
  double dist = 0.0;
  dist = sqrt( pow(this->x - othercar.x,2) + pow(this->y - othercar.y,2) );
  return dist;
}

// determine the lane othercar is driving in
lane PathFinder::in_lane(double d) {
  if (d < 0.0) {
    return OUTSIDE_LANE_LEFT;
  } else if ((d >= 0.0) and (d <= 4.0)) {
    return LEFT_LANE;
  } else if ((d >  4.0) and (d <= 8.0)) {
    return MIDDLE_LANE;
  } else if ((d >  8.0) and (d <= 12.0)) {
    return RIGHT_LANE;
  } else {
    return OUTSIDE_LANE_RIGHT;
  }
}

// test lane detection
/*
for (int k=0; k<sensor_fusion.size(); k++) {
  double my_d = sensor_fusion[k][6];
  lane my_lane = pf.in_lane(my_d);
  cout << k << " d = " << my_d << " lane = " << my_lane << endl;
}
*/




/*

std::function<int(double)> f(char);
https://stackoverflow.com/questions/31387238/c-function-returning-function

lesson 5.30

def to_equation(coefficients):
    """
    Takes the coefficients of a polynomial and creates a function of
    time from them.
    """
    def f(t):
        total = 0.0
        for i, c in enumerate(coefficients): 
            total += c * t ** i
        return total
        return f

*/



/*

def transition_function(predictions, current_fsm_state, current_pose, cost_functions, weights):
    # only consider states which can be reached from current FSM state.
    possible_successor_states = successor_states(current_fsm_state)

    # keep track of the total cost of each state.
    costs = []
    for state in possible_successor_states: 
        # generate a rough idea of what trajectory we would
        # follow IF we chose this state.
        trajectory_for_state = generate_trajectory(state, current_pose, predictions)

        # calculate the "cost" associated with that trajectory.
        cost_for_state = 0
        for i in range(len(cost_functions)) :
            # apply each cost function to the generated trajectory
            cost_function = cost_functions[i]
            cost_for_cost_function = cost_function(trajectory_for_state, predictions)

            # multiply the cost by the associated weight
            weight = weights[i]
            cost_for_state += weight * cost_for_cost_function
        costs.append({'state' : state, 'cost' : cost_for_state})

    # Find the minimum cost state.
    best_next_state = None
    min_cost = 9999999
    for i in range(len(possible_successor_states)):
        state = possible_successor_states[i]
        cost  = costs[i]
        if cost < min_cost:
            min_cost = cost
            best_next_state = state 

            return best_next_state

*/  
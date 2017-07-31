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

// pass a vector of vectors x to cout
void PathFinder::output_vector2(vector<vector<double>> x) {
  for (int i=0; i<x.size(); i++) {
    cout << "i=" << i << ": ";
    for (int j=0; j<x[j].size(); j++) {
          cout << x[i][j] << " ";
    }
    cout << endl; 
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

// determine index of closest othercar in lane (in front)
int PathFinder::distance2car_inlane(vector<Vehicle> othercars, double my_s, double my_d) {
  int othercar_id = -1;
  
  // smallest distance
  double sd = 99999;
  double coll_s  = -9999;
  lane coll_lane;
  
  // get my lane
  lane my_lane = in_lane(my_d);
  
  for (int i=0; i < othercars.size(); i++) {
    
    // get distance of car
    double distance = distance2car(othercars[i]);
    
    // get lane of car
    lane othercarlane = in_lane(othercars[i].d);
    
    // get s of car
    double s = othercars[i].s;
    
    if ((othercarlane == my_lane) && (distance < sd) && (s > my_s)) {
      othercar_id = othercars[i].id;
      sd = distance;
      coll_lane = othercarlane;
      coll_s = s;
    }
    
  }
  //cout << "dist2car:= \t" << othercar_id << "\t s = " << coll_s << "\t d = " << my_d <<"\t dist = "<< sd<< endl;
        
  // -1 if no car detected
  return othercar_id;
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

// calculate state of target vehicle using offset delta at time T 
/*
     target_vehicle - id of leading vehicle (int) which can be used to retrieve
       that vehicle from the "predictions" dictionary. This is the vehicle that 
       we are setting our trajectory relative to.

     delta - a length 6 array indicating the offset we are aiming for between us
       and the target_vehicle. So if at time 5 the target vehicle will be at 
       [100, 10, 0, 0, 0, 0] and delta is [-10, 0, 0, 4, 0, 0], then our goal 
       state for t = 5 will be [90, 10, 0, 4, 0, 0]. This would correspond to a 
       goal of "follow 10 meters behind and 4 meters to the right of target vehicle"

     T - the desired time at which we will be at the goal (relative to now as t=0)
*/
vector<double> PathFinder::predictions0(Vehicle targetcar, double T, vector<double> delta) {
  // get current state of target vehicle & predict state of targetvehicle at time T
  vector<double> target_state = targetcar.state_at(T);
  // store in temporary vectors
  vector<double> a = delta; //{100, 10, 0, 0, 0, 0};
  vector<double> b = target_state; //{-10, 0, 0, 4, 0, 0}; 
  // add vector delta to it
  transform (a.begin(), a.end(), b.begin(), b.begin(), plus<double>()); // result in b!!!
  return b;
}
// test predictions
//double T = 5;
//vector<double> delta = {10,0,0,4,0,0};
//vector<double> target_state = pf.predictions(mytargetcar, T, delta);
//pf.output_vector(target_state);


// generate a bunch of samples using gaussion noise
vector<vector<double>> PathFinder::perturb_goal0(vector<double> target_state, int n_samples) { 
  // generate gaussians
  default_random_engine gen;
  // creates a normal (Gaussian) distribution for s, d
  normal_distribution<double> dist_s1_init(target_state[0], SIGMA_S[0]);
  normal_distribution<double> dist_s2_init(target_state[1], SIGMA_S[1]);
  normal_distribution<double> dist_s3_init(target_state[2], SIGMA_S[2]);
  normal_distribution<double> dist_d1_init(target_state[3], SIGMA_D[0]);
  normal_distribution<double> dist_d2_init(target_state[4], SIGMA_D[1]);
  normal_distribution<double> dist_d3_init(target_state[5], SIGMA_D[2]);  
  // create vector for samples
  vector<vector<double>> all_goals;
  for (int i=0; i<n_samples; i++) {
    // generate elements of new pertubed state
    double s            = dist_s1_init(gen);
    double s_dot        = dist_s2_init(gen);
    double s_double_dot = dist_s2_init(gen);
    double d            = dist_d1_init(gen);
    double d_dot        = dist_d2_init(gen);
    double d_double_dot = dist_d2_init(gen);
    // store in new state vector
    all_goals.push_back({s, s_dot, s_double_dot, d, d_dot, d_double_dot});    
  }
  return all_goals;
}

// PTG part 0: main function. that calls PTG 1-2 und gives back a new path
vector<vector<double>> PathFinder::PTG_0_main(vector<Vehicle> othercars, double velocity, vector<double> start_state, double horizon) { 
  
  // vector that gives back the generated path
  vector<vector<double>> mypath;
  
  // vector to keep the input to pertubate goals
  vector<vector<double>> input_goals;
  
  // vector to keep all_goals which will be evaluated by their costs
  vector<vector<double>> all_goals; // s, v, a d, d_dot, d_double_dot
  
  // vector to keep the_goal for state which will be added to all_goals
  vector<double> the_goal; // s, v, a d, d_dot, d_double_dot
  
  // vector that keeps the delta between
  vector<double> delta; // TODO: define          
  
  // determine s & d path increments
  
          // define maximum car speed
          double max_car_speed = 0.9 * SPEED_LIMIT; // apply safety factor of 90%
          
          // increment along s in m
          double b = 0.90*10.0; //10.0; // m/s2
          double c = 0.90*10.0; //50.0; // m/s3
          double dt = 0.02;        // s          
          double x0 = 0; //car_s;
          double v0 = max_car_speed; // change to variable speed according to traffic conditions later
          
          double vel_set;
          if (velocity < max_car_speed) {
            vel_set = velocity;
          } else {
            vel_set = max_car_speed;
          }
          
          double ds = x0 + vel_set * dt + b * dt*dt/2 + c * dt*dt*dt/6;
          
  // determine possible succesor states based on actual state
  vector<states> possible_successor_states = successor_states(state); 
  
  // loop over possible succesor states
  for (int istate=0; istate < possible_successor_states.size(); istate++) {
    
    // output of current investigated state
    if (be_verbose) cout << "investigating state: " << possible_successor_states[istate] << endl;;
            
    // check for state & add new_goals for investigated state to all_goals        
    switch (possible_successor_states[istate]) {
      case KL: 
        cout << "KL   possible" << endl;    
        // check whether a vehicles is close
        // ....
        
        // case 1 "vehicle close": set the_goal and store in all_goals
        
        
        // case 2 "no vehicle close": set the_goal and store in all_goals
        the_goal = {start_state[0] + ds , velocity, 0, start_state[3], 0, 0};
        input_goals.push_back(the_goal);
        
      case PLCL:
        cout << "PLCL possible" << endl;
        
      case PLCR:
        cout << "PLCR possible" << endl;
        
      case LCL:
        cout << "LCR  possible" << endl;
        // set the_goal and store in all_goals
        the_goal = {start_state[0] + ds , velocity, 0, start_state[3]-4, 0, 0};
        input_goals.push_back(the_goal);
        
      case LCR:
        cout << "LCR  possible" << endl;
        // set the_goal and store in all_goals
        the_goal = {start_state[0] + ds , velocity, 0, start_state[3]+4, 0, 0};
        input_goals.push_back(the_goal);
        
      default:
        cerr << "PathFinder::PTG_0_main: state not defined -> " << possible_successor_states[istate] << endl;
        break;
    } 
    
    // PTG part 1b: perturbate input_goals into all_goals >> generate a bunch of salternative (pertubed) goals using gaussion noise based on target_states during T...T-4*dt
    int n_samples = 20;
    double T = 2; // predict at time T missing !!!!!
    for (int i=0; i < input_goals.size(); i++) {
      //vector<vector<double>> new_goals = perturb_goal0(input_goals[i], n_samples);
      vector<vector<double>> new_goals = PTG_1b_all_goals(input_goals[i], T, n_samples, 0.5); 
      all_goals.insert(all_goals.end(),new_goals.begin(),new_goals.end()); // all_goals.push_back(new_goals);
    } 
    
    // trajectory_for_state = generate_trajectory(state, current_pose, predictions) --> all_goals
    // PTG part 2: generate trajectories for all_goals
    //vector<double> current_state = {pf.s,pf.v,pf.a,pf.d,0,0 }; // check this d_dot, d_double_dot!!!!
    vector<vector<double>> trajectories = PTG_2_trajectories(all_goals, start_state); // coefficients of trajectories!!!
    cout << "*** " << trajectories.size() << " new trajectories generated ***" << endl;
    
    // calculate the "cost" associated with that trajectory.
    // TODO: implement!......
    // take coefficients of trajectories == trajectories AND target_states == all_goals as input for cost functions
    
  }          
  
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
  
  return mypath;
}

// PTG part 1a: generate a bunch of salternative (pertubed) goals using gaussion noise based on target_states during T...T-4*dt
vector<vector<double>> PathFinder::PTG_1a_all_goals(Vehicle targetcar, double T, vector<double> delta, int n_goals, double dt) { 
  
  // create vector for samples
  vector<vector<double>> all_goals;
  
  // while t <= T + 4 * timestep:
  double t = T - 4 * dt;
  
  while (t <= T + 4 * dt) {
    
    // get current state of target vehicle & predict state of targetvehicle at time t
    vector<double> a = delta; 
    vector<double> b = targetcar.state_at(t); 
    // add vector delta to it
    transform (a.begin(), a.end(), b.begin(), b.begin(), plus<double>()); // result in b!!!
    vector<double> target_state = b;
    
    // generate gaussians
    default_random_engine gen;
    // creates a normal (Gaussian) distribution for s, d
    normal_distribution<double> dist_s1_init(target_state[0], SIGMA_S[0]);
    normal_distribution<double> dist_s2_init(target_state[1], SIGMA_S[1]);
    normal_distribution<double> dist_s3_init(target_state[2], SIGMA_S[2]);
    normal_distribution<double> dist_d1_init(target_state[3], SIGMA_D[0]);
    normal_distribution<double> dist_d2_init(target_state[4], SIGMA_D[1]);
    normal_distribution<double> dist_d3_init(target_state[5], SIGMA_D[2]);  
    
    // generate goals
    for (int i=0; i<n_goals; i++) {
      // generate elements of new pertubed state
      double s            = dist_s1_init(gen);
      double s_dot        = dist_s2_init(gen);
      double s_double_dot = dist_s2_init(gen);
      double d            = dist_d1_init(gen);
      double d_dot        = dist_d2_init(gen);
      double d_double_dot = dist_d2_init(gen);
      // store in new state vector
      all_goals.push_back({s, s_dot, s_double_dot, d, d_dot, d_double_dot, t});    
    }
    
    // increment time
    t += dt;
  }
  
  return all_goals;
}
/*
          // PTG part 1
          // generate a bunch of alternative (pertubed) goals using gaussion noise based on target_states during T...T-4*dt
          int n_goals = 50;
          double timestep = 0.5;
          vector<vector<double>> all_goals = pf.PTG_1_all_goals(mytargetcar, T, delta, n_goals, timestep);         
          cout << "*** " << all_goals.size() << " new goals generated ***" << endl;
          //vector<vector<double>> all_goals = pf.perturb_goal0(target_state, n_samples);
*/
    
// PTG part 1b: generate a bunch of salternative (pertubed) goals using gaussion noise based on target_states during T...T-4*dt
vector<vector<double>> PathFinder::PTG_1b_all_goals(vector<double> target_state, double T, int n_goals, double dt) { 
  
  // create vector for samples
  vector<vector<double>> all_goals;
  
  // while t <= T + 4 * timestep:
  double t = T - 4 * dt;
  
  // predict state at time T: s movement
  double fut_s     = this->v * t + this->a * t * t / 2;
  double fut_v     = this->a * t;
  target_state[0] += fut_s;
  target_state[1] += fut_v;
  cout << "PathFinder::PTG_1b_all_goals: predicting s" << endl;
      
  while (t <= T + 4 * dt) {
    
    // generate gaussians
    default_random_engine gen;
    
    // creates a normal (Gaussian) distribution for s, d
    normal_distribution<double> dist_s1_init(target_state[0], SIGMA_S[0]);
    normal_distribution<double> dist_s2_init(target_state[1], SIGMA_S[1]);
    normal_distribution<double> dist_s3_init(target_state[2], SIGMA_S[2]);
    normal_distribution<double> dist_d1_init(target_state[3], SIGMA_D[0]);
    normal_distribution<double> dist_d2_init(target_state[4], SIGMA_D[1]);
    normal_distribution<double> dist_d3_init(target_state[5], SIGMA_D[2]);  
    
    // generate goals
    for (int i=0; i<n_goals; i++) {
      // generate elements of new pertubed state
      double s            = dist_s1_init(gen);
      double s_dot        = dist_s2_init(gen);
      double s_double_dot = dist_s2_init(gen);
      double d            = dist_d1_init(gen);
      double d_dot        = dist_d2_init(gen);
      double d_double_dot = dist_d2_init(gen);
      // store in new state vector
      all_goals.push_back({s, s_dot, s_double_dot, d, d_dot, d_double_dot, t});    
    }
    
    // increment time
    t += dt;
  }
  
  return all_goals;
}      

// PTG part 2: generate trajectories for all_goals
vector<vector<double>> PathFinder::PTG_2_trajectories(vector<vector<double>> all_goals, vector<double> current_state) { 
  
  // create vector for trajectories
  vector<vector<double>> trajectories;
  
  // split current state into s and d component
  vector<double> start_s(3); 
  vector<double> start_d(3);
  vector<double> goal_s(3);
  vector<double> goal_d(3);
  
  start_s[0] = current_state[0];          
  start_s[1] = current_state[1];          
  start_s[2] = current_state[2];
  start_d[0] = current_state[3];          
  start_d[1] = current_state[4];          
  start_d[2] = current_state[5];  
   
  // loop over all_goals
  for (int i = 0; i < all_goals.size(); i++) {
    
    // get  data from all_goals
    goal_s[0] = all_goals[i][0];          
    goal_s[1] = all_goals[i][1];          
    goal_s[2] = all_goals[i][2];
    goal_d[0] = all_goals[i][3];          
    goal_d[1] = all_goals[i][4];          
    goal_d[2] = all_goals[i][5];
    double t = all_goals[i][6];  
   
    // generate coefficients for s & d
    vector<double> s_coefficients = PathFinder::JMT(start_s, goal_s, t);
    vector<double> d_coefficients = PathFinder::JMT(start_d, goal_d, t);
      
    //store reults
    //trajectories.push_back({s_coefficients,d_coefficients,t}); 
  trajectories.push_back({s_coefficients[0],s_coefficients[1],s_coefficients[2],d_coefficients[0],d_coefficients[1],d_coefficients[2],t});
    
  }  
  
  return trajectories;
}
/*          
          // PTG part 2
          // generate trajectories for all_goals
          vector<double> current_state = {pf.s,pf.v,pf.a,pf.d,0,0 }; // check this d_dot, d_double_dot!!!!
          vector<vector<double>> trajectories = pf.PTG_2_trajectories(all_goals, current_state);
          cout << "*** " << trajectories.size() << " new trajectories generated ***" << endl;
*/ 

// generate predictions of given allcars over horizon
//vector<vector<double>> PathFinder::CARpredictions(vector<Vehicle> mycars, int horizon = 10) {
map<int, vector<vector<double>>> PathFinder::CARpredictions(vector<Vehicle> mycars, int horizon = 10) {
    // https://stackoverflow.com/questions/32679740/in-c-how-to-insert-key-and-value-in-empty-map-from-another-full-map
    // http://www.cplusplus.com/reference/map/map/insert/
    // https://stackoverflow.com/questions/6952486/recommended-wagy-to-insert-elements-into-map
    // myMap[ key ] = value;
    //assert( myMap.find( key )->second == value ); // post-condition
    // map.insert(std::pair<key_type, value_type>(key, value));
    
    //if (be_verbose) cout << "PathFinder::CARpredictions..." << endl;;
      
    // store predictions here
    map<int, vector<vector<double>>>  predictions;
    
    // time step
    double dt = 0.2; 

    for (int j = 0; j < mycars.size(); j++) {
      
      // check for this car
      Vehicle mycar = mycars[j];
      vector<vector<double>> storevec;
      if (be_verbose) cout << "car id: " << mycar.id << endl;;
      
      // loop over horizon
      for( int i = 0; i < horizon; i++)
      {
        vector<double> carstate = mycar.state_at(i*dt); // {s, v, this->a, d, d_dot, this->d_double_dot}
        vector<double> store = {carstate[0], carstate[3]};
        storevec.push_back(store);
        //assert(predictions.find( mycar.id )->second == store ); // crashed während runtime
        //if (be_verbose) cout << i << " " << carstate[0] << " " << carstate[3] << endl;
      }
      predictions.insert(pair<int, vector<vector<double>>>(mycar.id, storevec));
    }  

    return predictions;

}

// get id of car with car_id in vector of Vehicles
int PathFinder::car_id(vector<Vehicle> cars, int car_id) {
  for (int i = 0; i < cars.size(); i++) {
    Vehicle car = cars[i];
    if (car.id == car_id) {
      return i;
    }
  }
  return -1; // not found
}

/*
// summarize all costs of individual cost functions
double PathFinder::cost_summary(vector<double> &traj, vector<double> &goal) {
  
}
*/


// Penalizes trajectories that span a duration which is longer or shorter than the duration requested.
double PathFinder::cost4duration(double t, double T) {
  return logistic(float(abs(t-T)) / T);
}

// costs for total acceleration
double PathFinder::cost4s_total_acc(vector<double> traj_coeff, vector<double> target_state, double dt, int horizon) {
  
  // split trajectory_coeff in s & d
  vector<double> s_coeff = {traj_coeff[0], traj_coeff[1], traj_coeff[2]};
  vector<double> d_coeff = {traj_coeff[3], traj_coeff[4], traj_coeff[5]};
  
  // differentiate polynomals twice
  vector<double> s_dot_coeff = differentiate(s_coeff);
  vector<double> d_dot_coeff = differentiate(d_coeff);
  vector<double> s_dot2_coeff = differentiate(s_dot_coeff);
  vector<double> d_dot2_coeff = differentiate(d_dot_coeff);
  
  // evaluate polynomals over horizon
  for (int i = 0; i < horizon; i++) {
    double t = float(i) * dt;
    double s_dot2 = evaluate_polynomal(s_dot2_coeff, t);
    double d_dot2 = evaluate_polynomal(d_dot2_coeff, t);
    if ((s_dot2 > MAX_ACCEL) || (d_dot2 > MAX_ACCEL)) {
      return 1;
    } 
  }
  return 0;
}

// costs for total jerk
double PathFinder::cost4s_total_jerk(vector<double> traj_coeff, vector<double> target_state, double dt, int horizon) {
  
  // split trajectory_coeff in s & d
  vector<double> s_coeff = {traj_coeff[0], traj_coeff[1], traj_coeff[2]};
  vector<double> d_coeff = {traj_coeff[3], traj_coeff[4], traj_coeff[5]};
  
  // differentiate polynomals 3x
  vector<double> s_dot_coeff = differentiate(s_coeff);
  vector<double> d_dot_coeff = differentiate(d_coeff);
  vector<double> s_dot2_coeff = differentiate(s_dot_coeff);
  vector<double> d_dot2_coeff = differentiate(d_dot_coeff);
  vector<double> s_dot3_coeff = differentiate(s_dot2_coeff);
  vector<double> d_dot3_coeff = differentiate(d_dot2_coeff);
    
  // evaluate polynomals over horizon
  for (int i = 0; i < horizon; i++) {
    double t = float(i) * dt;
    double s_dot3 = evaluate_polynomal(s_dot3_coeff, t);
    double d_dot3 = evaluate_polynomal(d_dot3_coeff, t);
    if ((s_dot3 > MAX_JERK) || (d_dot3 > MAX_JERK)) {
      return 1;
    } 
  }
  return 0;
}

/*
// Penalizes trajectories whose s coordinate (and derivatives) differ from the goal.
double PathFinder::cost4s_diff(vector<double> trajectory, vector<double> target_state, double t, double T) {
  //vector<double> 
  //eval: 
  // evaluate a polynomal 5th order determined by coefficients at value x
  // double PathFinder::evaluate_polynomal(vector<double> coeffcients, double x)
  
  // calculates the derivative of a polynomial and returns the corresponding coefficients. ref: helpers.py in lesson 5.30
  // vector<double> PathFinder::differentiate(vector<double> coefficients)
  
  //double s = trajectory[0];
  
  return -1; //logistic(float(abs(t-T)) / T);
}
*/

/*
def s_diff_cost(traj, target_vehicle, delta, T, predictions):
    """
    Penalizes trajectories whose s coordinate (and derivatives) 
    differ from the goal.
    """
    s, _, T = traj
    target = predictions[target_vehicle].state_in(T)
    target = list(np.array(target) + np.array(delta))
    s_targ = target[:3]
    S = [f(T) for f in get_f_and_N_derivatives(s, 2)]
    cost = 0
    for actual, expected, sigma in zip(S, s_targ, SIGMA_S):
        diff = float(abs(actual-expected))
        cost += logistic(diff/sigma)
        return cost

*/
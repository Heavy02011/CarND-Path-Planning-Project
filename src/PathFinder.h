#ifndef PATHFINDER_H
#define PATHFINDER_H

#include <random>
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
#include <map>
//#include "Dense"

#include "vehicle.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

struct test_case {
	
		vector<double> start;
		vector<double> end;
		double T;
};

// data structure for FSM following the ideas of https://stackoverflow.com/questions/133214/is-there-a-typical-state-machine-implementation-pattern

typedef enum { KL,           // 0 keep lane
               PLCL,         // 1 prepare lane change left
               PLCR,         // 2 prepare lane change right
               LCL,          // 3 lane change left
               LCR } states; // 4 lane change right             
 
typedef enum { OUTSIDE_LANE_LEFT,
               LEFT_LANE,
               MIDDLE_LANE,
               RIGHT_LANE,
               OUTSIDE_LANE_RIGHT } lane;
               
class PathFinder 
{

  public:
        
    PathFinder(double x, double y, double d, double s, double v, double a, double yaw);

    virtual ~PathFinder();

    // test output
    bool be_verbose = true;
        
    // number of sample points on track
    const int N_SAMPLES = 50;
      
    // convert from MPH to m/s
    const double MPH2MPS = 0.44704;
    
    // set max speed im m/s
    const double SPEED_LIMIT = 50.0 * MPH2MPS;
    
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
    //VectorXd SIGMA_S = VectorXd(3);
    //VectorXd SIGMA_D = VectorXd(3);
    vector<double> SIGMA_S = { 10.0, 4.0, 2.0 }; // SIGMA for s, s_dot, s_double_dot
    vector<double> SIGMA_D = { 1.0, 1.0, 1.0 };  // SIGMA for d, d_dot, d_double_dot
    const double SIGMA_T = 2.0;
    
    // define initial state
    states state = KL;
    
    // define a vector of all states
    vector<states> all_states = { KL, PLCL, PLCR, LCL, LCR }; 
    
    //###
    
    // definition of lesson 4.16
  
    double d;
  
    double x;
    
    double y;
    
    double s;
  
    double v;
  
    double a;
    
    double yaw;
      
    double target_speed;
  
    //int lanes_available;
  
    double max_acceleration;
  
    double goal_lane;
  
    double goal_s;
    
    // id of target vehicle to follow
    int target_vehicle; // 99 if no vehicle in range
    
    
    //###
    
    // sample function
    double myfunc(double deg);
    
    // Calculate the Jerk Minimizing Trajectory that connects the initial state to the final state in time T.
    vector<double> JMT(vector< double> start, vector <double> end, double T);
    
    // generate straight path
    void generate_path_straight(double car_x, double car_y, double car_yaw, vector<double> *next_x, vector<double> *next_y); //, vector<double> *prev_x, vector<double> *prev_y);
 
    // generate straight path
    void generate_path_circle(double car_x, double car_y, double car_yaw, vector<double> *next_x, vector<double> *next_y, vector<double> *prev_x, vector<double> *prev_y);
    
    // define 
    vector<states> successor_states(states input);
    
    //===== helper functions =====
    
    // A function that returns a value between 0 and 1 for x in the 
    // range [0, infinity] and -1 to 1 for x in the range [-infinity, infinity].
    // Useful for cost functions.
    double logistic(double x);
    
    // pass a vector x to cout
    void output_vector(vector<double> x);
    
    // pass a vector of vectors x to cout
    void output_vector2(vector<vector<double>> x);
    
    // evaluate a polynomal 5th order determined by coefficients at value x
    double evaluate_polynomal(vector<double> coeffcients, double x);
    
    //std::function<int()> f();
    //std::function<int(double)> f(char);
    
    //function<double()> to_equation();
    //function<double(vector<double>)> to_equation(vector<double> coefficents);
    
    //calculates the derivative of a polynomial and returns the corresponding coefficients. ref: helpers.py in lesson 5.30
    vector<double> differentiate(vector<double> coefficients);
     
    // calculate distance of othercar to my car 
    double distance2car(Vehicle othercar);
    
    // calculate distance of othercar to my car (in front)
    int distance2car_inlane(vector<Vehicle> othercar, double my_s, double my_d);
    
    // determine the lane othercar is driving in
    lane in_lane(double d);
    
    // calculate state of target vehicle using offset delta at time T 
    vector<double> predictions0(Vehicle targetcar, double T, vector<double> delta);
    
    // generate a bunch of samples using gaussion noise
    vector<vector<double>> perturb_goal0(vector<double> target_state, int n_samples);
   
    // PTG part 0: main function. that calls PTG 1-2 und gives back a new path
    vector<vector<double>> PTG_0_main(vector<Vehicle> othercars, double velocity, vector<double> start_state, double horizon);
    
    // PTG part 1a: generate a bunch of salternative (pertubed) goals using gaussion noise based on target_states during T...T-4*dt
    vector<vector<double>> PTG_1a_all_goals(Vehicle targetcar, double T, vector<double> delta, int n_goals, double dt);
    
    // PTG part 1b: generate a bunch of salternative (pertubed) goals using gaussion noise based on target_states during T...T-4*dt
    vector<vector<double>> PTG_1b_all_goals(vector<double> target_state, double T, int n_goals, double dt);
    
    // PTG part 2: generate trajectories for all_goals
    vector<vector<double>> PTG_2_trajectories(vector<vector<double>> all_goals, vector<double> current_state);
    
    //vector<vector<double>> CARpredictions(vector<Vehicle>, int horizon);
    map<int, vector<vector<double>>> CARpredictions(vector<Vehicle>, int horizon);
    
    // get id of car with car_id in vector of Vehicles
    int car_id(vector<Vehicle> cars, int car_id);
    
    // summarize all costs of individual cost functions
    double cost_summary(vector<double> traj_coeff, vector<double> target_state, double dt, int horizon, double t, double T, vector<Vehicle> othercars);
    
    // Penalizes trajectories that span a duration which is longer or shorter than the duration requested.
    double cost4duration(double t, double T);
    
    // costs for total acceleration
    double cost4total_acc(vector<double> traj_coeff, vector<double> target_state, double dt, int horizon);
    
    // costs for total jerk
    double cost4total_jerk(vector<double> traj_coeff, vector<double> target_state, double dt, int horizon);
    
    // costs for deviation from d
    double cost4d_diff(vector<double> traj_coeff, vector<double> target_state, double dt, int horizon);
    
    // costs for driving not close to target velocity
    double cost4v_diff(vector<double> traj_coeff, vector<double> target_state, double dt, int horizon);
    
    // costs for driving not close to target velocity
    double cost4collision(vector<double> traj_coeff, vector<Vehicle> othercars, double dt, int horizon);
    
    // save mypath
    void savepath(string file_path, string vname, vector<vector<double>> mypath, int n_elements);
    
    vector<vector<double>> newpath(vector<vector<double>> myprevious_path_x, vector<vector<double>> myprevious_path_y, vector<vector<double>> mysensor_fusion);
  private:
      
      
};




#endif /* PATHFINDER_H */
#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h" // rbx
#include "PathFinder.h"

#include <cmath>
#include <vector>

#include "vehicle.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace tk; //rbx

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
  
//rbx
  
  // output waypoints & size 
  cout << "waypoints.size: " << map_waypoints_x.size() << std::endl;
  for (int i=0; i<map_waypoints_x.size(); i++) {
    cout << i <<  " map_waypoints_x,y,s,dx,dy = " << map_waypoints_x[i] << " ";  
    cout << map_waypoints_y[i] << " ";     
    cout << map_waypoints_s[i] << " ";  
    cout << map_waypoints_dx[i] << " "; 
    cout << map_waypoints_dy[i] << endl;
  }
  
  // setup spline objects
  spline waypointspline_x;
  spline waypointspline_y;
  spline waypointspline_dx;
  spline waypointspline_dy;
      
  // generate splines
  waypointspline_x.set_points( map_waypoints_s, map_waypoints_x );
  waypointspline_y.set_points( map_waypoints_s, map_waypoints_y );
  waypointspline_dx.set_points(map_waypoints_s, map_waypoints_dx);
  waypointspline_dy.set_points(map_waypoints_s, map_waypoints_dy);
  
  // get new point at my_s with --> double y = waypointspline_x(my_s);
  //double x_smooth = waypointspline_x(100.0);
  //cout <<  x_smooth << endl;
            

  // generate PathFinder class & initial states
  //PathFinder pf;
  PathFinder pf(0, 0, 0, 0, 0, 0, 0); // double x, double y, double d, double s, double v, double a, double yaw
  
  pf.state = KL;
  
  // setup a counter for every time step
  int counter = 0;
  
  /*
    h.onMessage([&count,&pp, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&WP_spline_x,&WP_spline_y,&WP_spline_dx,&WP_spline_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
    uWS::OpCode opCode) {
  */
  
//rbx
  //ur h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
      
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &waypointspline_x, &waypointspline_y, &waypointspline_dx, &waypointspline_dy, &pf, &counter](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
            
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
            // ["sensor_fusion"] A 2d vector of cars and then that car's 
            // [car's unique ID, car's x position in map coordinates, car's y position in map coordinates, car's x velocity in m/s, car's y velocity in m/s, car's s position in frenet coordinates, car's d position in frenet coordinates.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            
          /*                ["sensor_fusion"] A 2d vector of cars and then that car's 
          
                            car's unique ID, 
                            car's x position in map coordinates, 
                            car's y position in map coordinates, 
                            car's x velocity in m/s, 
                            car's y velocity in m/s, 
                            car's s position in frenet coordinates, 
                            car's d position in frenet coordinates.
                                                  
          sensor_fusion = [ [0 ,1022.598,1147.169,14.19477 ,10.3549   ,237.563  , 8.870655],
                            [1 ,1044.688,1155.066,15.69208 ,6.39629   ,261.0865 , 9.971673],
                            [2 ,1133.548,1187.77 ,14.78464 ,1.11656   ,357.6642 , 1.827328],
                            [3 ,904.5436,1124.697, 2.529505,2.679815  ,119.9534 ,10.10191 ],
                            [4 ,949.4061,1134.024,15.38097 ,1.409398  ,164.616  , 2.596563],
                            [5 ,1156.399,1189.374,15.09513 ,1.014567  ,380.5711 , 1.956437],
                            [6 ,1166.188,1181.94 ,16.5754  ,0.8770281 ,390.0564 , 9.954964],
                            [7 ,901.1926,1124.588, 2.646314,0.4060514 ,124.8682 ,10.17496 ],
                            [8 ,892.0728,1130.129,18.83123 ,1.271136  ,107.481  , 4.67399 ],
                            [9 ,1096.84 ,1174.888,17.36887 ,4.232598  ,332.6492 ,10.20708 ],
                            [10,1071.899,1174.909,16.43713 ,6.889896  ,293.7626 , 2.037921],
                            [11,836.1995,1132.793,20.70742 ,0.03320942, 51.60747, 2.123863] ]  
                      
          car x,y,s,d,yaw,speed = 909.48 1128.67 124.834 6.16483 0 0
                    
          */    
            
          
          
            
// =====================================================================================          
// ==== start implementation ===========================================================
// =====================================================================================
          
          cout << "***********************************************************" << endl;
          cout << "*** time step: " << counter << endl;
          cout << "*** car s,d,yaw,speed = " << car_s << " " << car_d << " " << car_yaw << " " << car_speed << endl;   
          cout << "***********************************************************" << endl;
          
          // output of path data  
          cout << "==========================================================================" << endl;          
          cout << "previous_path_x = " << previous_path_x << endl;  
          cout << "previous_path_y = " << previous_path_y << endl;     
          cout << "end_path_s = " << end_path_s << endl;  
          cout << "end_path_d = " << end_path_d << endl; 
          cout << "sensor_fusion = " << sensor_fusion << endl; 
          //cout << "*** car x,y,s,d,yaw,speed = " << car_x << " " << car_y << " " << car_s << " " << car_d << " " << car_yaw << " " << car_speed << endl; 
          cout << "*** car x,y = " << car_x << " " << car_y<< endl;         
          cout << endl;
          cout << "==========================================================================" << endl;          
          cout << endl;                
            
          // ***************************************************************************
          // 1 update my cars data in PathFinder object
          // ***************************************************************************
          
          // store actual car state & DELAY 1-2s ????????????????????  
          pf.s = car_s;
          pf.d = car_d;
          pf.x = car_x;
          pf.y = car_y;
          pf.v = car_speed*pf.MPH2MPS;
  

          // ***************************************************************************
          // 2 limit investigation of traffic to a smaller range out of sensor_fusion
          // ***************************************************************************
          
          // keep vehicles in range in this vector      
          vector<Vehicle> vehicles_inrange;
          vector<Vehicle> all_cars;
          
          // select the car that has largest distance & velocity in front of us as target
          
              // create instance for mytargetcar
              Vehicle mytargetcar(99,0,0,0,0,0,0,0,0,0);
              
              // largest distance in front
              double ld_front = -999999;
              
              // largest velocity in front
              double lv_front = -999999;
              
              // id of potential target car
              int id_front;
          
          // loop over all vehicles
          for (int k=0; k<sensor_fusion.size(); k++) {
            
            // generate a test vehicle for every entry in sensor_fusion
            Vehicle othercar(sensor_fusion[k][0], sensor_fusion[k][1], sensor_fusion[k][2], sensor_fusion[k][3], sensor_fusion[k][4], sensor_fusion[k][5], sensor_fusion[k][6], 0, 0, 0);
            
            // determine distance to my car
            double distance = pf.distance2car(othercar);
            //cout << "distance = " << distance << endl;
            
            // check whether car is suitable as target car withen next 50-500m in front of us
            if ((distance > ld_front) and ((othercar.s > car_s+50) and (othercar.s < car_s+500)) and (othercar.v > lv_front)) {
              id_front = k;
              mytargetcar = othercar;
            }
              
            // store vehicles in range
            if (distance < 50) {
              vehicles_inrange.push_back(othercar);
            }
            all_cars.push_back(othercar);
            
          }    
          // show vehicles_inrange 
          cout << "*** " << vehicles_inrange.size() << " vehicles in range 50m detected ***" << endl;
          for (int ii=0; ii<vehicles_inrange.size(); ii++) {
            vehicles_inrange[ii].display(vehicles_inrange[ii]);
          }
        
          // calculate predictions of vehicles_inrange over time horizon time steps
          int horizon = 50;
          map<int, vector<vector<double>>> predictions = pf.CARpredictions(vehicles_inrange, horizon);      
          cout << "*** " << predictions.size() << " predictions generated ***" << endl;
       
          // output of map predictions
          map<int, vector<vector<double>>>::iterator it= predictions.begin();
          while(it != predictions.end())
          {
              int car_id = it->first;
              vector<vector<double>> state_vector = it->second;
              cout << "key: " << car_id << " " << endl;
              pf.output_vector2(state_vector); 
              it++;
          }
          
            

          // ***************************************************************************
          // 3 select a real or virtual vehicle to follow & set target state
          // ***************************************************************************  
            
          //TODO: ...complete this... 
          cout << "*** id of target car: " << id_front << " ***" << endl;
          mytargetcar.display(mytargetcar);
                
          // predict state of target vehicle at time T using offset delta
          double T = 5;
          vector<double> delta = {10,0,0,4,0,0};
          vector<double> target_state = pf.predictions0(mytargetcar, T, delta);
          cout << "target_state" << endl;
          pf.output_vector(target_state);
          
          // PTG part 1
          // generate a bunch of alternative (pertubed) goals using gaussion noise based on target_states during T...T-4*dt
          int n_goals = 50;
          double timestep = 0.5;
          vector<vector<double>> all_goals = pf.PTG_1_all_goals(mytargetcar, T, delta, n_goals, timestep);         
          cout << "*** " << all_goals.size() << " new goals generated ***" << endl;
          //vector<vector<double>> all_goals = pf.perturb_goal0(target_state, n_samples);
          
          // PTG part 2
          // generate trajectories for all_goals
          vector<double> current_state = {pf.s,pf.v,pf.a,pf.d,0,0 }; // check this d_dot, d_double_dot!!!!
          vector<vector<double>> trajectories = pf.PTG_2_trajectories(all_goals, current_state);
          cout << "*** " << trajectories.size() << " new trajectories generated ***" << endl;
          
          
          
          
          // ***************************************************************************
          // 4 run finite state machine & evaluate possible trajectories & costs
          // ***************************************************************************

          //TODO: implement PTG.py here
          
          
/*
def transition_function(predictions, current_fsm_state, current_pose, cost_functions, weights):
    
          # only consider states which can be reached from current FSM state.
          possible_successor_states = successor_states(current_fsm_state)
*/
          
          // determine possible succesor states based on actual state
          vector<states> possible_successor_states = pf.successor_states(pf.state); 

/*
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
     
          // loop over possible succesor states
          for (int istate=0; istate < possible_successor_states.size(); istate++) {
            
            // output of current investigated state
            cout << "investigating state: " << possible_successor_states[istate] << endl;;
            
            // check for possible collision in the future and adjust costs
            lane my_lane = pf.in_lane(pf.d);
            //TODO: create a Vehicle for myself and make a prediction as for all_cars
            //....
            double costs = 0;
            map<int, vector<vector<double>>>::iterator it= predictions.begin();
            while(it != predictions.end())
            {
                int car_id = it->first;
                vector<vector<double>> state_vector = it->second; // vectors of s, d over horizon
                int othercar_id = pf.car_id(all_cars, car_id); // get index in all_cars for car_id
                double d = all_cars[othercar_id].d; // othercar with car_id
                double s = all_cars[othercar_id].s; // othercar s with car_id
                lane othercar_lane = pf.in_lane(d); // lane of the othercar
                
                //TODO:  function(state_vector) --> find collision index for s
                if ((othercar_lane == my_lane) and ((state_vector[49][0]-s)< 5)) {
                  cout << "/////////////////////////////////" << endl;
                  cout << "////// collision ahead //////////" << endl;
                  cout << "/////////////////////////////////" << endl; 
                  cout << "car_id: " << othercar_id << " s: " << s << " d: " << d <<endl;
                }
                
                //cout << "key: " << car_id << " " << endl;
                //pf.output_vector2(state_vector); 
                
                it++;
            }
            
            // check for state
            if (possible_successor_states[istate] == KL )  cout << "KL   possible" << endl;
            if (possible_successor_states[istate] == PLCL) cout << "PLCL possible" << endl;
            if (possible_successor_states[istate] == PLCR) cout << "PLCR possible" << endl;
            if (possible_successor_states[istate] == LCL)  cout << "LCR  possible" << endl;
            if (possible_successor_states[istate] == LCR)  cout << "LCR  possible" << endl;
            
            // generate path for current state
            // trajectory_for_state = generate_trajectory(state, current_pose, predictions)           
            vector<double> start = {car_s, 0, 0};
            vector<double> end   = {car_s, 0, 0};
            double T = 1;

            // setup input for JMT & output results
            start = {pf.s, pf.v, pf.a}; // actual state of my car
            end   = {pf.s+10, pf.SPEED_LIMIT, pf.MAX_ACCEL};
            T     = 20;
            //cout << "start/end/T = " << endl;
            //pf.output_vector(start);
            //pf.output_vector(end);
            //cout << T << endl;;
                        
            // call JMT & output results
            vector<double> coefficients = pf.JMT(start, end, T);
            //cout << "coefficients = " << endl;
            //pf.output_vector(coefficients);
            
            // test differentiate & output results
            vector<double> diff_coeff = pf.differentiate(coefficients);
            //cout << "diff_coefficients = " << endl;
            //pf.output_vector(diff_coeff);
            
            // determine costs according to lesson 4.16 behavioural planning for each successing possible state
            // #############################################################
            // ....
            
            
            
            
            
          }          
            
          // ***************************************************************************
          // 5 generate path coordinates
          // ***************************************************************************
          
          //TODO: ...complet this...
            
            
          // ***************************************************************************
          // XX generate path (test with circle)
          // ***************************************************************************
                  
          // variables for actual vehicle data
          double pos_x;
          double pos_y;
          double angle;
          
          // current size of path
          int path_size = previous_path_x.size();

          // store previous path in new one for smooth transition 
          for(int i = 0; i < path_size; i++)
          {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
          }

          // treat case of initial path of length zero, one and two
          if(path_size == 0)
          {
              pos_x = car_x;
              pos_y = car_y;
              angle = deg2rad(car_yaw);
          }
          else
          {
              pos_x = previous_path_x[path_size-1];
              pos_y = previous_path_y[path_size-1];

              double pos_x2 = previous_path_x[path_size-2];
              double pos_y2 = previous_path_y[path_size-2];
              angle = atan2(pos_y-pos_y2,pos_x-pos_x2);
          }
          
          next_x_vals.push_back(pos_x);
          next_y_vals.push_back(pos_y);
          // generate circle
          /*
          double dist_inc = 0.5;
          for(int i = 0; i < 50-path_size; i++)
          {    
              next_x_vals.push_back(pos_x+(dist_inc)*cos(angle+(i+1)*(pi()/100)));
              next_y_vals.push_back(pos_y+(dist_inc)*sin(angle+(i+1)*(pi()/100)));
              
              pos_x += (dist_inc)*cos(angle+(i+1)*(pi()/100));
              pos_y += (dist_inc)*sin(angle+(i+1)*(pi()/100));
            
          }
          */
          
          // ########
          
          //double max_car_speed = 0.8 * 50.0 * 0.44704; // apply safety factor of 90%
          double max_car_speed = 1.0 * pf.SPEED_LIMIT; // apply safety factor of 90%
          
          // position
          double b = 0.8*10.0; //10.0; // m/s2
          double c = 0.8*10.0; //50.0; // m/s3
          double dt = 0.02;        // s          
          double x0 = 0; //car_s;
          double v0 = max_car_speed; // change to variable speed according to traffic conditions later
          double xt = x0 + v0 * dt + b * dt*dt/2 + c * dt*dt*dt/6;
          cout << "==========================================================================" << endl;
          cout << "car_s, car_d, xt, v0 = " << car_s << " " << car_d << " " << xt << " " << v0 << endl;
          cout << endl;
          double x_smooth = waypointspline_x(car_s);
          cout <<  "smooth x = " << x_smooth << endl;
          cout << endl;
          cout << endl;
          cout << "==========================================================================" << endl;
          
          double ds = xt; //0.4; // increment along s for time interval dt  
          double dd = 2.0; //car_s - 6.0; //0.0; // increment for changing lane  
          
          // ########
          
          
          // generate path in middle lane
          double dist_inc = 0.5;
          //for(int i = 0; i < 50-path_size; i++)
          for(int i = path_size; i < 50; i++)              
          {  
           
            // new increments along s,d
            //double delta_s = car_s + i * ds;
            //double delta_d = car_d + i * dd;
            double delta_s = pf.s + i * ds;
            double delta_d = pf.d + i * dd;
            
            // use spline to get smooth new path points of road center
            double new_x0 = waypointspline_x(delta_s);
            double new_y0 = waypointspline_y(delta_s);
            double new_dx0 = waypointspline_dx(delta_s);
            double new_dy0 = waypointspline_dy(delta_s);
                       
            // adjust for lane
            double new_x = new_x0 + new_dx0 * dd;
            double new_y = new_y0 + new_dy0 * dd;
            
            // store the calculated path
            next_x_vals.push_back(new_x); 
            next_y_vals.push_back(new_y);   
           
           /* 
            // straight line
            double my_x = car_x+(dist_inc*i)*cos(deg2rad(car_yaw));
            double my_y = car_y+(dist_inc*i)*sin(deg2rad(car_yaw));
            
            next_x_vals.push_back(my_x); 
            next_y_vals.push_back(my_y);   
            */
          }
          
          
          // increment the counter
          counter += 1;
          
          
// =====================================================================================          
// ==== end implementation =============================================================
// =====================================================================================
          
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}


//########
/*
int NextWaypoint(double x, double y, double theta, 
                 vector<double> &maps_x, vector<double> &maps_y,
                 vector<double> &maps_dx, vector<double> &maps_dy)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

//	double heading = atan2( (map_y-y),(map_x-x) );
	
    //heading vector
    double hx = map_x-x;
    double hy = map_y-y;
    
    //Normal vector:
    double nx = maps_dx[closestWaypoint];
    double ny = maps_dy[closestWaypoint];
    
    //Vector into the direction of the road (perpendicular to the normal vector)
    double vx = -ny;
    double vy = nx;
    
    //Here we assume that the vehicle goes in the right directions
    //(If this is not the case, we have to examine theta and we might need to decrease closestWaypoint!)
    
    //If the inner product of v and h is positive then we are behind the waypoint so we do not need to
    //increment closestWaypoint, otherwise we are beyond the waypoint and we need to increment closestWaypoint.

    double inner = hx*vx+hy*vy;
    if (inner<0.0) {
        closestWaypoint++;
    }
    
//    double heading = atan2( vy,vx );

//	double angle = abs(theta-heading);

//	if(angle > pi()/4)
//	{
//		closestWaypoint++;
//	}
    

    return closestWaypoint;
}
*/
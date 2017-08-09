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
#include <string> 

#include <cmath>
#include <vector>

#include "vehicle.h"
//#include "matplotlibcpp.h"

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
/*  
  // output waypoints & size 
  cout << "waypoints.size: " << map_waypoints_x.size() << std::endl;
  for (int i=0; i<map_waypoints_x.size(); i++) {
    cout << i <<  " map_waypoints_x,y,s,dx,dy = " << map_waypoints_x[i] << " ";  
    cout << map_waypoints_y[i] << " ";     
    cout << map_waypoints_s[i] << " ";  
    cout << map_waypoints_dx[i] << " "; 
    cout << map_waypoints_dy[i] << endl;
  }
*/  
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
  
  // upsample waypoints ref: https://discussions.udacity.com/t/latency-handling/322156
  vector<double> map_waypoints_x_upsampled;
  vector<double> map_waypoints_y_upsampled;
  vector<double> map_waypoints_s_upsampled;
  //vector<double> map_waypoints_dx_upsampled;
  //vector<double> map_waypoints_dy_upsampled;
    
  spline spline_x, spline_y;
  spline_x.set_points(map_waypoints_s, map_waypoints_x);
  spline_y.set_points(map_waypoints_s, map_waypoints_y);

  // refine path with spline.
  int spline_samples = 12000;
  for (size_t i = 0; i < spline_samples; ++i) {
    map_waypoints_x_upsampled.push_back(spline_x(i));
    map_waypoints_y_upsampled.push_back(spline_y(i));
    map_waypoints_s_upsampled.push_back(i);
  }
  
  
  // get new point at my_s with --> double y = waypointspline_x(my_s);
  //double x_smooth = waypointspline_x(100.0);
  //cout <<  x_smooth << endl;
            

  // generate PathFinder class & initial states
  //PathFinder pf;
  PathFinder pf(0, 0, 0, 0, 0, 0, 0); // double x, double y, double d, double s, double v, double a, double yaw
  
  pf.state = KL;
  
  // setup a counter for every time step
  int counter = 0;

  // initialize velocity          
  double vel_set = 0; //49.5 * pf.MPH2MPS;

  // set my lane
  int my_lane = 1;


  
  /*
    h.onMessage([&count,&pp, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&WP_spline_x,&WP_spline_y,&WP_spline_dx,&WP_spline_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
    uWS::OpCode opCode) {
  */
  
//rbx
  //ur h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
      
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &waypointspline_x, &waypointspline_y, &waypointspline_dx, &waypointspline_dy, &pf, &counter, &map_waypoints_x_upsampled, &map_waypoints_y_upsampled, &map_waypoints_s_upsampled, &vel_set, &my_lane](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
            //vector<vector<double>> myprevious_path_x = j[1]["previous_path_x"]; 
            //vector<vector<double>> myprevious_path_y = j[1]["previous_path_y"];

          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
            // ["sensor_fusion"] A 2d vector of cars and then that car's 
            // [car's unique ID, car's x position in map coordinates, car's y position in map coordinates, car's x velocity in m/s, car's y velocity in m/s, car's s position in frenet coordinates, car's d position in frenet coordinates.
          	auto sensor_fusion = j[1]["sensor_fusion"];
            //vector<vector<double>> mysensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	//vector<double> next_x_vals;
          	//vector<double> next_y_vals;


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
                      
          cout << endl;             
          cout << "***********************************************************" << endl;
          cout << "*** time step: " << counter << endl;
          cout << "*** car s,d,yaw,speed = " << car_s << " " << car_d << " " << car_yaw << " " << car_speed << endl;   
          cout << "*** car_x, y = " << car_x << " " << car_y << endl;   
          cout << "***********************************************************" << endl;
          
          // ***************************************************************************
          // 0 new implementation of path
          // ***************************************************************************

          // get size of previous path
          int previous_path_size = previous_path_x.size();

          // make smooth transition
         
          if (previous_path_size > 0) {
            car_s = end_path_s;
          }

          // are other cars too close
          bool othercars_too_close = false;

          // determine speed to drive
          for (int i=0; i<sensor_fusion.size(); i++) {
            // car in my lane
            float d = sensor_fusion[i][6];
            if(d < (2+4*my_lane+2) && (2+4*my_lane-2)) {
              double vx = sensor_fusion[i][3];
              double vy = sensor_fusion[i][4];
              double check_speed = sqrt(vx*vx+vy*vy);
              double check_car_s = sensor_fusion[i][5];

              // predict cars future state
              check_car_s += previous_path_size*0.02*check_speed;
              if ((check_car_s > car_s) && ((check_car_s - car_s) < 30)) {
                //vel_set = 29.5*pf.MPH2MPS;
                othercars_too_close = true;

                // just change lane left
                //if (my_lane > 0) {
                //  my_lane = 0;
                //}

              }

            }
          }

          // adjust velocity
          if (othercars_too_close) {
            vel_set -= 0.224*pf.MPH2MPS;
          } else if (vel_set < 49.5*pf.MPH2MPS) {
            vel_set += 0.224*pf.MPH2MPS;
          }


          // ***************************************************************************
          // 1 update my cars data in PathFinder object
          // ***************************************************************************

          pf.s = car_s;
          pf.d = car_d;
          pf.x = car_x;
          pf.y = car_y;
          pf.v = car_speed*pf.MPH2MPS;
          pf.yaw = car_yaw;
          double my_dx = waypointspline_dx(car_s); //TODO: upsampling necessary???????????????
          double my_dy = waypointspline_dy(car_s);

         
          // place my car in a vector of cars as well to be able to use state prediction
          vector<Vehicle> mycars;
          Vehicle mycar(23, car_x, car_y, my_dx*car_speed*pf.MPH2MPS, my_dy*car_speed*pf.MPH2MPS, car_s, car_d, 0, 0, 0);
          mycars.push_back(mycar); // this will be done every time step again
  
          // predict the state of my car in the future (1s)
          vector<double> mystate = mycar.state_at(1);
          vector<double> car_xy = getXY(mystate[0], mystate[3], map_waypoints_s, map_waypoints_x, map_waypoints_y);
/*          
          // upddate car data to future state
          car_s = mystate[0];
          car_d = mystate[3];
          car_x = car_xy[0];
          car_y = car_xy[1];
*/          
          // {s, v, this->a, d, d_dot, this->d_double_dot};

          // save all relevant input data
          // vector<vector<double>> mynewpath = pf.newpath(myprevious_path_x, myprevious_path_y, mysensor_fusion);

          // ***************************************************************************
          // 2 limit investigation of traffic to a smaller range out of sensor_fusion
          // ***************************************************************************
          
          // keep vehicles in range in this vector      
          vector<Vehicle> cars_inrange;
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
              cars_inrange.push_back(othercar);
            }
            all_cars.push_back(othercar);
            
          }    
          // show cars_inrange 
          cout << "*** " << cars_inrange.size() << " vehicles in range 50m detected ***" << endl;
          for (int ii=0; ii<cars_inrange.size(); ii++) {
              cout << "cars_inrange = \t" << cars_inrange[ii].id << "\t s = " << all_cars[cars_inrange[ii].id].s <<"\t d = " << all_cars[cars_inrange[ii].id].d << "\t dist = " << pf.distance2car(all_cars[cars_inrange[ii].id]) << endl;
                      
              cars_inrange[ii];
          }           
/*          
          if (pf.be_verbose) {
            cout << "*** " << cars_inrange.size() << " vehicles in range 50m detected ***" << endl;
            for (int ii=0; ii<cars_inrange.size(); ii++) {
              cars_inrange[ii].display(cars_inrange[ii]);
            }
          }
*/                      

          // ***************************************************************************
          // 3 select a real or virtual vehicle to follow & set target state
          // ***************************************************************************  
            
          //TODO: ...complete this... 
/*
          if (pf.be_verbose) {
            cout << "*** id of target car: " << id_front << " ***" << endl;
            mytargetcar.display(mytargetcar);
          }
                
          // predict state of target vehicle at time T using offset delta
          double T = 5;
          vector<double> delta = {10,0,0,4,0,0};
          vector<double> target_state = pf.predictions0(mytargetcar, T, delta);
          
          // output of target state
          if (pf.be_verbose) {
            cout << "target_state" << endl;
            pf.output_vector(target_state);
          }
*/          
          
          // ***************************************************************************
          // 4 run finite state machine & evaluate possible trajectories & costs
          // ***************************************************************************
/*
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
*/          
  
          // predict the state of my car in the future (1s)
          //mystate = mycar.state_at(horizon*0.02);
          mystate = mycar.state_at(2);
                    
          // check for possible collision in the future and adjust costs
          //lane my_lane = pf.in_lane(mystate[3]); //pf.in_lane(pf.d);   // s, v, this->a, d, d_dot, this->d_double_dot
//double my_s = car_s; //mystate[0]; #######################
//double my_d = car_d; //mystate[3];
double my_s = mystate[0]; 
double my_d = mystate[3];
          cout << "mystate: \t" << "\t s = " << car_s <<"\t d = " <<car_d<< endl;
          
          //get id of closest in front car in my lane 
          int closecar_id = pf.distance2car_inlane(all_cars, my_s, my_d); // TODO: limit to cars_inrange!!!!!!!
          
          // get distance of that car
          double closecar_dist = pf.distance2car(all_cars[closecar_id]);
          cout << "closecar = \t" << closecar_id << "\t s = " << all_cars[closecar_id].s <<"\t d = " << all_cars[closecar_id].d << "\t dist = " << closecar_dist << endl;
          cout << "myfutstate: \t" << "\t s = " << mystate[0] <<"\t d = " <<mystate[3]<< endl;
                    
          cout << endl;

/*          
          // crude test drive
          double pos_d = car_d; //6;
          int closecar_id_right  = pf.distance2car_inlane(all_cars, my_s, 10);
          int closecar_id_middle = pf.distance2car_inlane(all_cars, my_s, 6);
          int closecar_id_left   = pf.distance2car_inlane(all_cars, my_s, 2);
          double closecar_dist_right  = pf.distance2car(all_cars[closecar_id_right]);
          double closecar_dist_middle = pf.distance2car(all_cars[closecar_id_middle]);
          double closecar_dist_left   = pf.distance2car(all_cars[closecar_id_left]);
          
          if ((closecar_dist_left < 25) && (closecar_dist_middle < 25) && (closecar_dist_right < 25)) {
            // slow down
            //pf.SPEED_LIMIT *= 0.9;
            cout << "<<< SLOW DOWN >>>" << endl;
          } else {
            // change lane
            if ((closecar_dist < 25) && (closecar_dist_right > 25)) {
            pos_d = 10;
            } else if ((closecar_dist < 25) && (closecar_dist_left > 25)) {
              pos_d = 2;
            } else if ((closecar_dist < 25) && (closecar_dist_middle > 25)){
              pos_d = 6;
            } else {
              // keep lane
              lane my_lane = pf.in_lane(mystate[3]);
              if (my_lane == LEFT_LANE) pos_d = 2; 
              if (my_lane == MIDDLE_LANE) pos_d = 6; 
              if (my_lane == RIGHT_LANE) pos_d = 10;           
            }
          }
                    
 */


          // ***************************************************************************
          // 5 generate path coordinates
          // ***************************************************************************         
                  
                    
          // Option 3: generate new path using walkthrough video approach

          // create widely spaced path points to construct the final path
          vector<double> ptsx;
          vector<double> ptsy;

          // define the cars reference state: starting point or previous path end point
          double ref_x = car_x;
          double ref_y = car_y;
          double ref_yaw = deg2rad(car_yaw);

          // if previous path size is only 1 point or smaller use car as starting reference
          cout << "previous_path_size = " << previous_path_size << endl;
          if (previous_path_size < 2) {
            // use two points to align car path with yaw
            double prev_car_x = car_x - cos(ref_yaw); //???????????????
            double prev_car_y = car_y - sin(ref_yaw); // ref_yaw???????
            ptsx.push_back(prev_car_x);
            ptsx.push_back(car_x);
            ptsy.push_back(prev_car_y);
            ptsy.push_back(car_y);

          } else
          // use last previous point as starting reference 
          {
            // set last point of previous path as reference state
            ref_x = previous_path_x[previous_path_size-1];
            ref_y = previous_path_y[previous_path_size-1];

            double ref_x_prev = previous_path_x[previous_path_size-2];
            double ref_y_prev = previous_path_y[previous_path_size-2];
            ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

            // save the two points
            ptsx.push_back(ref_x_prev);
            ptsx.push_back(ref_x);
            ptsy.push_back(ref_y_prev);
            ptsy.push_back(ref_y);

          }

          // get my lane (its in int my_lane)
          //lane my_lane = pf.in_lane(mystate[3]);
          cout << "my_lane = " << my_lane << endl;
          //int my_lane = 1;

          // generate the wide spaced ANCHOR points to generate a spline for the path
          vector<double> next_wp0 = getXY(car_s+30,(2+4*my_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp1 = getXY(car_s+60,(2+4*my_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp2 = getXY(car_s+90,(2+4*my_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

          // add ANCHOR points to vector
          ptsx.push_back(next_wp0[0]);
          ptsx.push_back(next_wp1[0]);
          ptsx.push_back(next_wp2[0]);
          ptsy.push_back(next_wp0[1]);
          ptsy.push_back(next_wp1[1]);
          ptsy.push_back(next_wp2[1]);

          // convert ANCHOR points for spline to local car coordinate system
          for (int i=0; i<ptsx.size(); i++) {
            // convert to local car coordinates and rotate yaw to zero degrees
            double shift_x = ptsx[i] - ref_x;
            double shift_y = ptsy[i] - ref_y;

            ptsx[i] = shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw);
            ptsy[i] = shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw);
            cout << i << " " << ptsx[i] << " " << ptsy[i] << endl;
          }

          // construct spline
          tk::spline s;
          s.set_points(ptsx, ptsy);

          // define the actual path points
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          // part 1: add previous path points
          for (int i = 0; i < previous_path_size; i++)
          {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          // break ANCHOR spline into segments. distance = N * 0.02 * velocity
          double target_x = 30.0; // set target 30m ahead
          double target_y = s(target_x);
          double target_dist = sqrt(target_x*target_x + target_y*target_y);

          // spline constructor x for getting y in local car coordinate system
          double x_add_on = 0;
          double N = target_dist / (0.02*vel_set);  // check for units!!!!!!!!!!!
          double dx = target_x/N;

          // part 2: add the new points by using the spline
          for (int i = 1; i < 50-previous_path_size; ++i)
          {
            // generate points in local coordinate system
            double x_point = x_add_on + dx;
            double y_point = s(x_point);

            // set new x
            x_add_on = x_point;

            // transfrom points back to global coordinate system
            double x_delta = x_point;
            double y_delta = y_point;
            x_point = ref_x + x_delta * cos(ref_yaw) - y_delta * sin(ref_yaw);
            y_point = ref_y + x_delta * sin(ref_yaw) + y_delta * cos(ref_yaw);

            // store points in vector
            next_x_vals.push_back(x_point);
            next_y_vals.push_back(y_point);

          }        
         
          // ***************************************************************************
          

         
          
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
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
  PathFinder pf;
  pf.state = KL;
  
  /*
    h.onMessage([&count,&pp, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&WP_spline_x,&WP_spline_y,&WP_spline_dx,&WP_spline_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
    uWS::OpCode opCode) {
  */
  
//rbx
  //ur h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
      
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &waypointspline_x, &waypointspline_y, &waypointspline_dx, &waypointspline_dy, &pf](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
            
// ==== rbx =======================================================================
                    
          /*
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
          */    
            
          // output of path data  
          cout << "==========================================================================" << endl;          
          cout << "previous_path_x = " << previous_path_x << endl;  
          cout << "previous_path_y = " << previous_path_y << endl;     
          cout << "end_path_s = " << end_path_s << endl;  
          cout << "end_path_d = " << end_path_d << endl; 
          cout << "sensor_fusion = " << sensor_fusion << endl;
          cout << "car x,y,s,d,yaw,speed = " << car_x << " " << car_y << " " << car_s << " " << car_d << " " << car_yaw << " " << car_speed << endl;       
          cout << endl;
          cout << "==========================================================================" << endl;          
          cout << endl;
          
          
          // ***************************************************************************
          // run finite state machine
          // ***************************************************************************
          
          cout << "test Finite State Machine" << endl;  
          
          // loop over all states
          for (int istate=0; istate<pf.all_states.size(); istate++) {
            cout << istate << " " << pf.all_states[istate] << endl;
            
            // set a new state
            pf.state = pf.all_states[istate];
            
            // check for possible succesor states
            vector<states> my_possible_states = pf.successor_states(pf.state); 
            
            // print results
            for (int kk=0; kk<my_possible_states.size(); kk++) {
              cout << my_possible_states[kk] << endl;
            } 
            
            // test check for state
            if (pf.state == LCL) cout << "LCR detected" << endl;
            
          }  
          
          
          // ***************************************************************************
          // generate path (test with circle)
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

          // generate circle
          double dist_inc = 0.5;
          for(int i = 0; i < 50-path_size; i++)
          {    
              next_x_vals.push_back(pos_x+(dist_inc)*cos(angle+(i+1)*(pi()/100)));
              next_y_vals.push_back(pos_y+(dist_inc)*sin(angle+(i+1)*(pi()/100)));
              
              pos_x += (dist_inc)*cos(angle+(i+1)*(pi()/100));
              pos_y += (dist_inc)*sin(angle+(i+1)*(pi()/100));
            
          }
          
          
          
          
// ==== rbx =======================================================================
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
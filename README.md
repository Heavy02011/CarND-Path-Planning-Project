# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program
   
## Rubic Points

### A The car is able to drive at least 4.32 miles without incident.
![<Display Name>](https://github.com/Heavy02011/P1-CarND-term3-Path-Planning-Project/blob/master/png/onelap.png)

### B There is a reflection on how to generate paths.
  
The following steps are taken to generate the path:

##### 1 update my cars data in PathFinder object

An object PathFinder is created to hold all relevant data of the car, the sensor fusion data and helping functions like predicting the future state of cars, getting distances between our car and traffic etc.

##### 2 limit investigation of traffic to a smaller range out of sensor_fusion

To save computational time and keep the path planner responsive the investigations are limited to traffic within a range of 50m.

##### 3 get id's and distance of close cars
This part of the path planner takes the latency of the simulator into account. Our cars state is predicted at a future state that is calculated as follows:

    double timelag = (50-previous_path_size)*0.02;
    
50 points are generated in total for the path. Every cycle only some points are 'consumed' by the simulator giving back a number of previous_path_size points from the previous cycle. The difference times the time increment of 0.02s gives the lagging time. The state of our car is then predicted at lagging time in the future.

The id's and distance of the cars around us are calculated for every lane like

    // cars with their id and distance around us
    int closecar_id_right  = pf.distance2car_inlane(all_cars, my_s, 10);
    int closecar_id_middle = pf.distance2car_inlane(all_cars, my_s, 6);
    int closecar_id_left   = pf.distance2car_inlane(all_cars, my_s, 2);
    double closecar_dist_right  = pf.distance2car(all_cars[closecar_id_right]);
    double closecar_dist_middle = pf.distance2car(all_cars[closecar_id_middle]);
    double closecar_dist_left   = pf.distance2car(all_cars[closecar_id_left]);

##### 4 run finite state machine & evaluate possible trajectories (simple model)
This part mainly follows the approach of the [Udacity Q&A walkthrough] (https://classroom.udacity.com/nanodegrees/nd013/parts/6047fe34-d93c-4f50-8336-b70ef10cb4b2/modules/27800789-bc8e-4adc-afe0-ec781e82ceae/lessons/23add5c6-7004-47ad-b169-49a5d7b1c1cb/concepts/3bdfeb8c-8dd6-49a7-9d08-beff6703792d): "Man, that was a monster!" This was the most difficult part of the project because the behavior of the simulator was not that clear.

The basic idea is to generate the path in a local coordinate system (x pointing along the s and y perpendicular to the right) and just use a few ANCHOR points to generate a spline with them. The starting points are taken from the previous path points preserving the heading of the car and at stations of 30, 60 and 90 m ahead.

The state machine uses the following states & decisions:
* set a bool othercars_too_close to true if cars are closer than 30m at the future state,
* a save lane change is regarded as possible if the distance to other cars is not smaller than 10m and there is enough distance (40m) on the neighbor lanes,
* an escape condition is set if we are stuck in heavy traffic,
* in all instances of heavy traffic around us the reduction of speed has priority,
* reduction of speed is limited to a minimum value of 30 mph.

A full JMT and PTG with cost functions was implemented in the PathFinder but currently not used here.

An alternative approach would have been a much simpler set of cost functions ranking the 'best' lane like described in this [article](https://medium.com/@mohankarthik/path-planning-in-highways-for-an-autonomous-vehicle-242b91e6387d).

However, after spending so much time on the project I preferred to stay with this simpler approach that just worked.

##### 5 generate path coordinates
Well, this follows exactly the approach as showed in [Udacity Q&A walkthrough] (https://classroom.udacity.com/nanodegrees/nd013/parts/6047fe34-d93c-4f50-8336-b70ef10cb4b2/modules/27800789-bc8e-4adc-afe0-ec781e82ceae/lessons/23add5c6-7004-47ad-b169-49a5d7b1c1cb/concepts/3bdfeb8c-8dd6-49a7-9d08-beff6703792d).

The basic idea is to generate the new points in a local car coordinate system and add them to the list of previous points. The total number of points is aimed at a total number of 50. After transforming them back to global map coordinates they are handed over to the simulator.


### C Final Remark
This was by far the most time consuming project of all terms, but could have been done in much less time as the walkthrough showed us. However, implementing a full PTG with JMT was a good practice. A test may follow after the submission.

----
### Simulator. You can download the Term3 Simulator BETA which contains the Path Planning Project from the [releases tab](https://github.com/udacity/self-driving-car-sim/releases).

In this project your goal is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. You will be provided the car's localization and sensor fusion data, there is also a sparse map list of waypoints around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 50 m/s^3.

#### The map of the highway is in data/highway_map.txt
Each waypoint in the list contains  [x,y,s,dx,dy] values. x and y are the waypoint's map coordinate position, the s value is the distance along the road to get to that waypoint in meters, the dx and dy values define the unit normal vector pointing outward of the highway loop.

The highway's waypoints loop around so the frenet s value, distance along the road, goes from 0 to 6945.554.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./path_planning`.

Here is the data provided from the Simulator to the C++ Program

#### Main car's localization Data (No Noise)

["x"] The car's x position in map coordinates

["y"] The car's y position in map coordinates

["s"] The car's s position in frenet coordinates

["d"] The car's d position in frenet coordinates

["yaw"] The car's yaw angle in the map

["speed"] The car's speed in MPH

#### Previous path data given to the Planner

//Note: Return the previous list but with processed points removed, can be a nice tool to show how far along
the path has processed since last time. 

["previous_path_x"] The previous list of x points previously given to the simulator

["previous_path_y"] The previous list of y points previously given to the simulator

#### Previous path's end s and d values 

["end_path_s"] The previous list's last point's frenet s value

["end_path_d"] The previous list's last point's frenet d value

#### Sensor Fusion Data, a list of all other car's attributes on the same side of the road. (No Noise)

["sensor_fusion"] A 2d vector of cars and then that car's [car's unique ID, car's x position in map coordinates, car's y position in map coordinates, car's x velocity in m/s, car's y velocity in m/s, car's s position in frenet coordinates, car's d position in frenet coordinates. 

## Details

1. The car uses a perfect controller and will visit every (x,y) point it recieves in the list every .02 seconds. The units for the (x,y) points are in meters and the spacing of the points determines the speed of the car. The vector going from a point to the next point in the list dictates the angle of the car. Acceleration both in the tangential and normal directions is measured along with the jerk, the rate of change of total Acceleration. The (x,y) point paths that the planner recieves should not have a total acceleration that goes over 10 m/s^2, also the jerk should not go over 50 m/s^3. (NOTE: As this is BETA, these requirements might change. Also currently jerk is over a .02 second interval, it would probably be better to average total acceleration over 1 second and measure jerk from that.

2. There will be some latency between the simulator running and the path planner returning a path, with optimized code usually its not very long maybe just 1-3 time steps. During this delay the simulator will continue using points that it was last given, because of this its a good idea to store the last points you have used so you can have a smooth transition. previous_path_x, and previous_path_y can be helpful for this transition since they show the last points given to the simulator controller with the processed points already removed. You would either return a path that extends this previous path or make sure to create a new path that has a smooth transition with this last path.

## Tips

A really helpful resource for doing this project and creating smooth trajectories was using http://kluge.in-chemnitz.de/opensource/spline/, the spline function is in a single hearder file is really easy to use.

---

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please (do your best to) stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html).

## Project Instructions and Rubric

Note: regardless of the changes you make, your project must be buildable using
cmake and make!


## Call for IDE Profiles Pull Requests

Help your fellow students!

We decided to create Makefiles with cmake to keep this project as platform
agnostic as possible. Similarly, we omitted IDE profiles in order to ensure
that students don't feel pressured to use one IDE or another.

However! I'd love to help people get up and running with their IDEs of choice.
If you've created a profile for an IDE that you think other students would
appreciate, we'd love to have you add the requisite profile files and
instructions to ide_profiles/. For example if you wanted to add a VS Code
profile, you'd add:

* /ide_profiles/vscode/.vscode
* /ide_profiles/vscode/README.md

The README should explain what the profile does, how to take advantage of it,
and how to install it.

Frankly, I've never been involved in a project with multiple IDE profiles
before. I believe the best way to handle this would be to keep them out of the
repo root to avoid clutter. My expectation is that most profiles will include
instructions to copy files to a new location to get picked up by the IDE, but
that's just a guess.

One last note here: regardless of the IDE used, every submitted project must
still be compilable with cmake and make./

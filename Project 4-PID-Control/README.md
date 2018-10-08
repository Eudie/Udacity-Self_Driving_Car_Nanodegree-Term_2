# Project: PID Controller
Self-Driving Car Engineer Nanodegree Program


---

### Overview
This is the 9th project of the [Self Driving Car Engineer Nanodegree](https://www.udacity.com/course/self-driving-car-engineer-nanodegree--nd013) I am taking part. <br>
The aim of this project was to build a PID controller in order to drive a car on a simulator. 


### Roles of P, I, and D
- **Proportional coefficient (P)**<br/>
This component had the most directly observable effect on the car's behavior. It causes the car to steer proportional to the car's distance from the lane center. If the car is far to the right it steers hard to the left, if it's slightly to the left it steers slightly to the right.
- **Integral coefficient (I)**<br/>
The role of this term is to account for systematic bias. If the vehicle has some manufacture error in steering then by integrating all error with time we can get sense of how our controller is doing. If it increases than this term will turn more to get on centre line. 
- **Differential coefficient (D)**<br/>
If we use only proportional coefficient than trajectory merge quickly to centre line but than overshoots. To make the transition smoother from CTE to centre line we use differential coefficient.


### Final Parameter

The final hyperparameters are P = 0.2, I = 0.0, D = 3.0, which I have selected using hit and trial in this version. Following videos are some example of my iterations:

[p=0.2,d=3.0,i=0-final_video](https://youtu.be/wbWURMcayuY)<br/>
[p=1,d=0,i=0](https://youtu.be/oBmG8EV3pAg)<br/>
[p=1,d=1,i=0](https://youtu.be/Q2rtHsRq2fc)

---

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1(mac, linux), 3.81(Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `./install-mac.sh` or `./install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```
    Some function signatures have changed in v0.14.x. See [this PR](https://github.com/udacity/CarND-MPC-Project/pull/3) for more details.
* Simulator. You can download these from the [project intro page](https://github.com/udacity/self-driving-car-sim/releases) in the classroom.

There's an experimental patch for windows in this [PR](https://github.com/udacity/CarND-PID-Control-Project/pull/3)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./pid`. 


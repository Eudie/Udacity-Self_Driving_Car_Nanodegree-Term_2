#define _USE_MATH_DEFINES

#include <cmath>

#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */

  x_ = F_*x_;
  //cout << F_ << endl;
	P_ = F_*P_*F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  // new state

  long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	MatrixXd S = H_*P_*H_.transpose() + R_;
	MatrixXd K = P_*H_.transpose()*S.inverse();


	x_ = x_ + K*(z-H_*x_);
	P_ = (I-K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  double rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  double theta = atan2(x_(1) , x_(0));
  double rho_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / rho;
  VectorXd h = VectorXd(3); // h(x_)
  h << rho, theta, rho_dot;

  long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);


	MatrixXd S = H_*P_*H_.transpose() + R_;
	MatrixXd K = P_*H_.transpose()*S.inverse();

  VectorXd y = z-h;


  if(y[1]>M_PI){
    // cout << y[1] <<endl;
    y[1]=y[1] - 2*M_PI;
    // cout << y[1] <<endl;
  }else if(y[1]< -1*M_PI){
    // cout << y[1] <<endl;
    y[1]=y[1] + 2*M_PI;
    // cout << y[1] <<endl;
  }
	x_ = x_ + K*(y);

	P_ = (I-K*H_)*P_;


}

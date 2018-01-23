#include "kalman_filter.h"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

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
  /*
    * predict the state
  */
	x_ = F_*x_;
	P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
	//std::cout<<"z size:"<<z.size()<<std::endl;
	VectorXd y;
	MatrixXd S;
	MatrixXd K;
	
	y = z - H_*x_; //error vec
	
	S = H_*P_*H_.transpose() + R_;
	K = P_*H_.transpose()*S.inverse();
	x_ += K*y;
	P_ = (Eigen::MatrixXd::Identity(x_.size(),x_.size()) - K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
	//convert x_ to measurement space
	//for clarity use some helper variables
	VectorXd y;
	MatrixXd S;
	MatrixXd K;
	
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);

	float dist_p = sqrt(px*px + py*py);

	VectorXd z_pred(3);
	if (dist_p > 1E-20)
		z_pred << dist_p, atan2(py, px), (px*vx + py*vy) / dist_p;
	else
		z_pred << dist_p, atan2(py, px), 0;
	y = z - z_pred; //error vec
	
	float y1_wrapped = fmod(y(1)+M_PI, 2*M_PI);
	if (y1_wrapped<0)
		y1_wrapped += 2*M_PI;
	y1_wrapped -= M_PI;
	//cout<<"zwrapped: "<<z_pred1_wrapped<<endl;
	//cout<<"z: "<<z_pred[1]<<endl;
	y[1] = y1_wrapped;

	S = H_*P_*H_.transpose()+ R_;
	if (fabs(S.determinant())>1E-20)// a better check would be to look at condition number
		K = P_*H_.transpose()*S.inverse();
	else
		K = MatrixXd::Zero(4,3);
	
	x_ += K*y;
	P_ = (MatrixXd::Identity(x_.size(), x_.size()) - K*H_)*P_;  
 }


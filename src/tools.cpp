#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /*
    * Calculate the RMSE here.
  */
	VectorXd ret_v = VectorXd::Zero(4);
	//error checking

	if (estimations.size() != ground_truth.size() || estimations.size() == 0) 
	{
		std::cout << "Error in input size" << std::endl;
		return ret_v;
	}
	for (int i = 0; i<estimations.size(); i++)
	{
		ret_v += VectorXd((estimations[i] - ground_truth[i]).array()*(estimations[i] - ground_truth[i]).array());
	}
	return ret_v;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /*
    * Calculate a Jacobian here.
  */
	MatrixXd ret_H = MatrixXd::Zero(3,4);
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//div by zero check
	if (px == 0 && py == 0)
	{
		std::cout << "div by zero" << std::endl;
		return ret_H;
	}

	//compute the Jacobian matrix
	float px2_py2 = px*px + py*py;
	ret_H(0, 0) = px / sqrt(px2_py2);
	ret_H(0, 1) = py / sqrt(px2_py2);
	ret_H(1, 0) = -py / (px2_py2);
	ret_H(1, 1) = px / (px2_py2);
	ret_H(2, 0) = py*(vx*py - vy*px) / pow(px2_py2, 1.5);
	ret_H(2, 1) = px*(vy*px - vx*py) / pow(px2_py2, 1.5);
	ret_H(2, 2) = px / sqrt(px2_py2);
	ret_H(2, 3) = py / sqrt(px2_py2);

	return ret_H;
}

#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.01, 0,
              0, 0.01;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0004, 0,
        0, 0, 0.09;
  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_ = tools.CalculateJacobian(VectorXd::Zero(4));

  ekf_.P_ =  Eigen::MatrixXd::Identity(4, 4);

  //measurement covariance laser
  ekf_.R_ = R_laser_;
  //measurement matrix
  ekf_.H_ = H_laser_;

  //the initial transition matrix F_
  ekf_.F_ = Eigen::MatrixXd::Identity(4, 4);
  //set the acceleration noise components
  noise_ax = 9;
  noise_ay = 9;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
    */
    ekf_.F_ = Eigen::MatrixXd::Identity(4, 4);  
   // first measurement
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	//measurements for Radar is (rho,phi,rho_dot)
	VectorXd xvec = VectorXd(measurement_pack.raw_measurements_);
	
	//ekf_.x_<< xvec[0] * cos(xvec[1]), xvec[0] * sin(xvec[1]), xvec[2] * cos(xvec[1]), xvec[2] * sin(xvec[1]);		
	ekf_.x_<< xvec[0] * cos(xvec[1]), xvec[0] * sin(xvec[1]), 0,0;		
	Hj_ = tools.CalculateJacobian(ekf_.x_);
	ekf_.P_ = MatrixXd::Identity(4,4)- Hj_.transpose()*((MatrixXd::Identity(3,3)+R_radar_).inverse())*Hj_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	VectorXd xvec = VectorXd(measurement_pack.raw_measurements_);
	ekf_.x_ << xvec[0], xvec[1], 0, 0;
	ekf_.P_ =  MatrixXd::Identity(4,4)-H_laser_.transpose()*((H_laser_*(H_laser_.transpose())+R_laser_).inverse())*H_laser_;

    }

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  //The STM F
  float del_t = (measurement_pack.timestamp_ - previous_timestamp_)*1E-6;
//  cout<<"delta_t_0: "<<del_t<<endl;
//  cout<<"time stamp: "<<previous_timestamp_<<"endl";
  if (del_t > 1)//probably started a new dataset, reset ic
  {
	del_t = 0;
	is_initialized_ = false;
  }
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_(0, 2) = del_t;
  ekf_.F_(1, 3) = del_t;

  float del_t2 = del_t*del_t;
  float del_t4 = del_t2*del_t2;
  float del_t3 = del_t2*del_t;  

//  cout<<"delta_t:"<<del_t<<endl;
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << del_t4 / 4 * noise_ax, 0, del_t3 / 2 * noise_ax, 0,
	  0, del_t4 / 4 * noise_ay, 0, del_t3 / 2 * noise_ay,
	  del_t3 / 2 * noise_ax, 0, del_t2*noise_ax, 0,
	  0, del_t3 / 2 * noise_ay, 0, del_t2*noise_ay;

  ekf_.Predict();
 
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	ekf_.R_ = R_radar_;
	ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	
  } else {

    // Laser updates
 	ekf_.R_ = R_laser_;
     	ekf_.H_ = H_laser_;
	ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}


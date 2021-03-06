#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
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
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

	VectorXd x = VectorXd(4);
	x << 0,0,0,0;

	MatrixXd P = MatrixXd(4,4);
	P << 	10, 0, 0, 0,
				0, 10, 0, 0,
				0, 0, 10, 0, 
				0, 0, 0, 10;
	
	MatrixXd F = MatrixXd(4,4);
	F <<  1, 0, 1, 0,
				0, 1, 0, 1,
				0, 0, 1, 0,
				0, 0, 0, 1;

	H_laser_ << 1,0,0,0,
							0,1,0,0;

	MatrixXd Q = MatrixXd(4,4);
	Q << 0, 0, 0, 0,
			 0, 0, 0, 0,
			 0, 0, 0, 0,
			 0, 0, 0, 0;
	
	// Initialize kalman filter with P and F
	// x H R and Q are dummy here, we dont know which sensor data comes first
	ekf_.Init(x,P,F,H_laser_,R_laser_,Q);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    ekf_.x_ = VectorXd(4);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
			double rho = measurement_pack.raw_measurements_[0];
			double phi = measurement_pack.raw_measurements_[1];
			double rho_dot = measurement_pack.raw_measurements_[2];
			// 
			ekf_.x_ << rho * cos(phi) , rho * sin(phi) , rho_dot * cos(phi) , rho_dot * sin(phi); 
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
			// 
			ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1],0,0;
    }
		previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = measurement_pack.timestamp_;

	ekf_.F_(0,2) = dt;
	ekf_.F_(1,3) = dt;

	double noise_ax = 9;
	double noise_ay = 9;

	ekf_.Q_ = MatrixXd(4,4);
	ekf_.Q_ << (noise_ax * pow(dt,4))/4, 0,	(noise_ax * pow(dt,3))/2,0,
							0, (noise_ay * pow(dt,4))/4,0,(noise_ay * pow(dt,3))/2,
							(noise_ax * pow(dt,3))/2, 0, noise_ax * pow(dt,2), 0,
							0, (noise_ay * pow(dt,3))/2, 0, noise_ay * pow(dt,2);

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
		ekf_.H_ = tools.CalculateJacobian(ekf_.x_);	
		ekf_.R_ = R_radar_;
		ekf_.UpdateEKF(measurement_pack.raw_measurements_);
		 
  } else {
    // TODO: Laser updates
		ekf_.H_ = H_laser_;
		ekf_.R_ = R_laser_;
		ekf_.Update(measurement_pack.raw_measurements_);
  }
}

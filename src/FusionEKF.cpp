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
   * Setting the process and measurement noises
   */
  
  // Initializing Laser Measurement Matrix 
  H_laser_ << 1, 0, 0, 0,
	     0, 1, 0, 0;

  // Initializing State Transition Matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
             
  // Initializing State Covariance Matrix
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

  // Initializing noise values for the Process Covariance matrix
  noise_ax = 9;
  noise_ay = 9;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessFirstMeasurement(const MeasurementPackage &measurement_pack) {
  /**
     * Initializing the state ekf_.x_ with the first measurement.
     */

  // first measurement
  cout << "EKF: " << endl;
  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 0.8, 0.6, 5.1, 1.8;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Convert radar from polar to cartesian coordinates 
    //         and initialize state. I am referring to phi as theta.
    float rho = measurement_pack.raw_measurements_(0);
    float theta = measurement_pack.raw_measurements_(1);
    float rhodot = measurement_pack.raw_measurements_(2);

    // Initialize state
    ekf_.x_(0) = rho * cos(theta);
    ekf_.x_(1) = rho * sin(theta);
    ekf_.x_(2) = rhodot * cos(theta);
    ekf_.x_(3) = rhodot * sin(theta);

    cout << "Radar First Measurement completed" << endl;
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    // Initialize state.
    ekf_.x_(0) = measurement_pack.raw_measurements_(0);
    ekf_.x_(1) = measurement_pack.raw_measurements_(1);
    cout << "Laser First Measurement Completed" << endl;
  }

  previous_timestamp_ = measurement_pack.timestamp_;

  // done initializing, no need to predict or update
  is_initialized_ = true;
}

void FusionEKF::InitializeEkfForPrediction(const MeasurementPackage &measurement_pack) {
  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Pre-computed values
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Initialize the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
         0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
         dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
         0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
}

void FusionEKF::InitializeEkfForUpdate(const MeasurementPackage &measurement_pack) {
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
  }
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    ProcessFirstMeasurement(measurement_pack);
    return;
  }

  /**
   * Prediction
   */
  InitializeEkfForPrediction(measurement_pack);
  ekf_.Predict();
  cout << "Prediction completed" << endl;

  /**
   * Update
   */
  InitializeEkfForUpdate(measurement_pack);
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_.Update(measurement_pack.raw_measurements_);
  }
  cout << "Update completed" << endl;

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

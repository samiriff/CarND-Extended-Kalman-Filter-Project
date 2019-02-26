#ifndef FusionEKF_H_
#define FusionEKF_H_

#include <fstream>
#include <string>
#include <vector>
#include "Eigen/Dense"
#include "kalman_filter.h"
#include "measurement_package.h"
#include "tools.h"

class FusionEKF {
 public:
  /**
   * Constructor.
   */
  FusionEKF();

  /**
   * Destructor.
   */
  virtual ~FusionEKF();

  /**
   * Run the whole flow of the Kalman Filter from here.
   */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
   * Kalman Filter update and prediction math lives in here.
   */
  KalmanFilter ekf_;

 private:
  /**
   * Processes the very first measurement and initializes state and timestamp
   */
  void ProcessFirstMeasurement(const MeasurementPackage &measurement_pack);

  /**
   * Initializes the State transition and Process Covariance matrices which will be used by the Kalman filter during prediction
   */
  void InitializeEkfForPrediction(const MeasurementPackage &measurement_pack);

  /**
   * Initializes the Measurement and Measurement Covariance matrices which will be used by the Kalman filter during updates
   */
  void InitializeEkfForUpdate(const MeasurementPackage &measurement_pack);

  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  Eigen::MatrixXd H_laser_;
  Eigen::MatrixXd Hj_;

  float noise_ax;
  float noise_ay;
};

#endif // FusionEKF_H_

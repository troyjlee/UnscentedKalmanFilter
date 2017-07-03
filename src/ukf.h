#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(const MeasurementPackage &meas_package);

    ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  // variable to hold the radar normalized innovation squared
  double NIS_radar_;

  // variable to hold the laser normalized innovation squared
  double NIS_laser_;

private:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  // lidar measurement matrix
  Eigen::Matrix<double,2,5> H_;

  // time of previous measurement
  long long previous_timestamp_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Covariance matrix of laser measurements 
  Eigen::Matrix<double,2,2> R_las_;

  // Covariance matrix of radar measurements 
  Eigen::Matrix<double,3,3> R_rad_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  // VectorXd weights_;
  // Try fixed size initialization
  Eigen::Matrix<double,15,1> weights_;

  ///* State dimension
  const int n_x_ = 5;

  ///* Augmented state dimension
  int n_aug_;

  // Number of sigma points
  int n_sigma_;

  ///* weight on mean sigma point
  double W0_;

  // Sigma scaling factor
  double sigma_scale_;

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double dt);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const VectorXd &z);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const VectorXd &z);

  // function for normalizing angles into [-pi,pi]
  void normalizeAngle(double &angle);

};

#endif /* UKF_H */

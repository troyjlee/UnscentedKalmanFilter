#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // No measurement read yet
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  
  //n_x_ = 5;
  n_aug_ = 7;
  n_sigma_ = 2*n_aug_ + 1;
  
  // Weight on the mean and sigma points scaling factor.
  W0_ = -1.0/3.0;
  sigma_scale_ = sqrt(n_aug_/(1-W0_));

  // initial state vector
  x_ = VectorXd(n_x_);

  // matrix to hold sigma point predictions
  Xsig_pred_ = MatrixXd(n_x_,n_sigma_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  // look into different initialization here
  P_ = MatrixXd::Identity(n_x_,n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .5;

  // Laser measurement covariance matrix
  R_las_ << 0.15*0.15, 0,
            0, 0.15*0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  R_rad_ = MatrixXd::Zero(3,3);
  R_rad_(0,0) = std_radr_ * std_radr_;
  R_rad_(1,1) = std_radphi_ * std_radphi_;
  R_rad_(2,2) = std_radrd_ * std_radrd_;
  
  H_ = MatrixXd::Zero(2,n_x_);
  H_(0,0) = 1;
  H_(1,1) = 1;

  // definition of weights to calculate
  // mean and covariance from predicted sigma points
  weights_(0) = W0_;
  for(int i=1; i < n_sigma_;++i){
      weights_(i) = (1-W0_)/(2*n_aug_);
  }
}

UKF::~UKF() {}

// helper function to keep angles in [-pi,pi]
void UKF::normalizeAngle(double &angle){
  while(angle > M_PI){
    angle -= 2*M_PI;
  }
  while(angle < -M_PI){
    angle += 2*M_PI;
  }
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  // Initialization step
  if (!is_initialized_) {
    cout << "UKF: " << endl;
    // initialize x vector
    // whether lidar or radar, velocity initialized to zero
    x_ = VectorXd::Zero(n_x_);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];
      x_(0) = rho*cos(phi);
      x_(1) = rho*sin(phi);
      x_(3) = phi;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = measurement_pack.raw_measurements_[0];
      x_(1) = measurement_pack.raw_measurements_[1];
    }


  previous_timestamp_ = measurement_pack.timestamp_;
  // we are now initialized
  is_initialized_ = true;
  return;
  }

  // compute the time elapsed between the current and previous measurements
  long long current_time = measurement_pack.timestamp_;
  // dt - expressed in seconds
  double dt = (current_time - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = current_time;
  
  // do prediction
  Prediction(dt);
  // do update

 if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_){
    UpdateLidar(measurement_pack.raw_measurements_);
  } 
 else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
    UpdateRadar(measurement_pack.raw_measurements_);
  } 
  cout << "x =" << endl;
  cout << x_ << endl;
  cout << "P_ =" << endl;
  cout << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */

void UKF::Prediction(double dt){
  // Generate Sigma points
  MatrixXd Xsig = MatrixXd(n_x_, n_sigma_);
  MatrixXd Xaug = MatrixXd(n_aug_, n_sigma_);
  Xsig = x_.replicate(1,n_sigma_);  

  // Cholesky factorization of covariance matrix P
  MatrixXd A = P_.llt().matrixL();
  for(int i = 0;i< n_x_;++i){
    Xsig.col(i+1) += sigma_scale_*A.col(i);
    Xsig.col(n_aug_+i+1) -= sigma_scale_*A.col(i);
  }
  // Add on block corresponding to noise variables
  MatrixXd Block = MatrixXd::Zero(2,n_sigma_);
  Block(0,n_x_+1) = sigma_scale_*std_a_;
  Block(1,n_x_+2) = sigma_scale_*std_yawdd_;
  Block(0,n_sigma_-2) = -sigma_scale_*std_a_;
  Block(1,n_sigma_-1) = -sigma_scale_*std_yawdd_;
  
  // vstack Xsig on top of Block
  Xaug << Xsig, Block;

  // Now we have the sigma points, predict them forward.
  for(int i=0;i < n_sigma_;++i){
    VectorXd u = Xaug.col(i);
    double v = u(2);
    double psi = u(3);
    double psi_dot = u(4);
    double noise_a = u(5);
    double noise_yaw = u(6);
    double dt2 = dt*dt;
    VectorXd noise = VectorXd(5);
    noise << noise_a*cos(psi)*dt2/2, noise_a*sin(psi)*dt2/2, dt*noise_a,
             dt2*noise_yaw/2, dt*noise_yaw;
    Xsig_pred_.col(i) = u.topRows(5) + noise;
    if(fabs(psi_dot) < 0.00001){
        Xsig_pred_(0,i) += v*cos(psi)*dt; 
        Xsig_pred_(1,i) += v*sin(psi)*dt; 
    } // endif
    else{
        Xsig_pred_(0,i) += v*(sin(psi+psi_dot*dt)-sin(psi))/psi_dot; 
        Xsig_pred_(1,i) += v*(-cos(psi+psi_dot*dt)+cos(psi))/psi_dot; 
        Xsig_pred_(3,i) += psi_dot*dt;
    } // endelse
  } // endfor

  // From predicted sigma points, generate mean and covariance
  // state weighted average
  //cout << "predicted sigma points" << endl;
  //cout << Xsig_pred_.row(3) << endl;
  x_ = Xsig_pred_ * weights_;
  
  // subtract the mean
  MatrixXd Xsig_meanzero(n_x_,n_sigma_);
  Xsig_meanzero = Xsig_pred_.colwise() - x_;
  for(int i=0;i < n_sigma_;++i){
    normalizeAngle(Xsig_meanzero(3,i));
  }
  // predict state covariance matrix
  P_ = Xsig_meanzero*weights_.asDiagonal()*Xsig_meanzero.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const VectorXd &z) {
  MatrixXd Ht = H_.transpose();
  MatrixXd S = MatrixXd(2,2);
  S = H_*P_*Ht + R_las_;
  MatrixXd K = MatrixXd(n_x_,2);
  K = P_*Ht*S.inverse();
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K*H_)*P_;
  VectorXd y(2);
  y = z - H_*x_;
  x_ = x_ + K*y;
  normalizeAngle(x_(3));
  NIS_laser_= y.transpose()*S.inverse()*y;
  ofstream lasfile("NIS_laser.txt",ios::app);
  lasfile << NIS_laser_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 */

void UKF::UpdateRadar(const VectorXd &z) {
  // first we compute measurement prediction from sigma points
  MatrixXd Zsig(3,n_sigma_);
  for(int i =0;i < n_sigma_; ++i){
      double px = Xsig_pred_(0,i);
      double py = Xsig_pred_(1,i);
      double v = Xsig_pred_(2,i);
      double psi = Xsig_pred_(3,i);
      double r = sqrt(px*px + py*py);
      if(r < 0.0000001){
        cout << "Error in UpdateRadar: Divide by zero" << endl;
      }
      else{
        Zsig(0,i) = r;
        Zsig(1,i) = atan2(py,px);
        Zsig(2,i) = (px*v*cos(psi) + py*v*sin(psi))/r;
      }
    }
    Eigen::Matrix<double,3,1> z_pred;
    z_pred = Zsig*weights_;
    MatrixXd Zsig_meanzero(3,n_sigma_);
    Zsig_meanzero = Zsig.colwise() - z_pred;
    for(int i=0;i < n_sigma_;++i){
      normalizeAngle(Zsig_meanzero(1,i));
    }
    MatrixXd S(3,3);
    S = Zsig_meanzero * weights_.asDiagonal() * Zsig_meanzero.transpose();
    S += R_rad_;

    // now we update based on predicted and actual measurement
    //calculate cross correlation matrix
    MatrixXd Xsig_meanzero(n_x_,n_sigma_);
    Xsig_meanzero = Xsig_pred_.colwise() - x_;
    for(int i=0;i < n_sigma_;++i){
      normalizeAngle(Xsig_meanzero(3,i));
    }
    MatrixXd Tc(n_x_,3);
    Tc = Xsig_meanzero * weights_.asDiagonal() * Zsig_meanzero.transpose();
    //calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    //update state mean and covariance matrix
    VectorXd y = z - z_pred;
    normalizeAngle(y(1));
    x_ += K*y;
    normalizeAngle(x_(3));
    P_ -= K*S*K.transpose();
    NIS_radar_ = y.transpose()*S.inverse()*y;
    ofstream radfile("NIS_radar.txt",ios::app);
    radfile << NIS_radar_ << endl;
}



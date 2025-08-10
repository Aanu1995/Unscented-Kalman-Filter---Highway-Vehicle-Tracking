#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // State dimension
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */


  time_us_ = 0.0;

  // Augmented state dimension
  n_aug_ = n_x_ + 2;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  // set weights
  weights_(0) = lambda_/(lambda_ + n_aug_);
  weights_.segment(1, 2 * n_aug_).setConstant(0.5 / (lambda_ + n_aug_));

  // Lidar Measurement noise covariance matrix
  R_lidar_ = Eigen::MatrixXd(2, 2);
  R_lidar_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  // Radar Measurement noise covariance matrix
  R_radar_ = Eigen::MatrixXd(3, 3);
  R_radar_ <<  std_radr_ * std_radr_, 0, 0,
        0, std_radphi_ * std_radphi_, 0,
        0, 0, std_radrd_ * std_radrd_;

  // The current NIS for lidar
  NIS_lidar_ = 0.0;

  // The current NIS for radar
  NIS_radar_ = 0.0;

  is_initialized_ = false;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * Switch between lidar and radar measurements.
   */

   // Initialize measurements on first call
    if(!is_initialized_) {

      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        double rho = meas_package.raw_measurements_(0); // range
        double phi = meas_package.raw_measurements_(1); // angle
        double rho_dot = meas_package.raw_measurements_(2); // radial velocity

        double cal_x = rho * cos(phi);
        double cal_y = rho * sin(phi);

        double cal_vx = rho_dot * cos(phi);
        double cal_vy = rho_dot * sin(phi);

        double cal_v = sqrt(cal_vx * cal_vx + cal_vy * cal_vy);

        x_ << cal_x , cal_y, cal_v, 0.0, 0.0;

        P_ << std_radr_ * std_radr_, 0, 0, 0, 0,
              0, std_radr_ * std_radr_, 0, 0, 0,
              0, 0, std_radrd_ * std_radrd_, 0, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1;

      } else {
        x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0.0, 0.0, 0.0;
        P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
              0, std_laspy_ * std_laspy_, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1;
      }

      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;

     return;
   }


  // calculate change in time between the current and previous update
   double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
   time_us_ = meas_package.timestamp_;

   // Prediction step
   Prediction(dt);

   // Update step
   if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
     UpdateLidar(meas_package);
   } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
     UpdateRadar(meas_package);
   }
}

void UKF::Prediction(double delta_t) {
  /**
   * Estimate the object's location. Modify the state vector, x_. 
   * Predict sigma points, the state, and the state covariance matrix.
   */

  // Predict the sigma points
  PredictSigmaPoints(delta_t);

  // Predict the state vector and state covariance matrix
  PredictStateMeanAndCovariance();
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  // set lidar measurement dimension
  int n_z = 2;

  // set true measurement from Lidar
  Eigen::VectorXd z = Eigen::VectorXd(n_z);
  z << meas_package.raw_measurements_(0),
       meas_package.raw_measurements_(1);

  // create matrix for sigma points in measurement space
  Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  Eigen::VectorXd z_pred = Eigen::VectorXd(n_z);
  z_pred.fill(0.0);

  // measurement covariance matrix S
  Eigen::MatrixXd S = Eigen::MatrixXd(n_z, n_z);
  S.fill(0.0);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    Zsig(0, i) =  Xsig_pred_(0, i);     // px
    Zsig(1, i) = Xsig_pred_(1, i);      // py
  }

  // mean predicted measurement
  z_pred.noalias() = Zsig * weights_;

  // innovation covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // adding noise
  S += R_lidar_;

  // create matrix for cross correlation Tc
  Eigen::MatrixXd Tc = Eigen::MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // normalize angles
    while (x_diff(3) > M_PI) x_diff(3) -= 2 * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  Eigen::MatrixXd K = Tc * S.inverse();

  // residual
  Eigen::VectorXd ZRealdiff = z - z_pred;

  // update state mean and covariance matrix
  x_ += K * ZRealdiff;
  P_ -= K * S * K.transpose();

  // Calculate Lidar NIS
  NIS_lidar_ = ZRealdiff.transpose() * S.inverse() * ZRealdiff;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  // set radar measurement dimension
  int n_z = 3;

  // set true measurement from Radar
  Eigen::VectorXd z = Eigen::VectorXd(n_z);
  z << meas_package.raw_measurements_(0),
       meas_package.raw_measurements_(1),
       meas_package.raw_measurements_(2);

   // create matrix for sigma points in measurement space
  Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  Eigen::VectorXd z_pred = Eigen::VectorXd(n_z);
  z_pred.fill(0.0);

  // measurement covariance matrix S
  Eigen::MatrixXd S = Eigen::MatrixXd(n_z, n_z);
  S.fill(0.0);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    // extract values for better readability
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v  = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);                          // r
    Zsig(1, i) = atan2(p_y, p_x);                                      // phi
    Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);  // r_dot
  }

  // mean predicted measurement
  z_pred.noalias() = Zsig * weights_;

  // innovation covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // adding noise
  S += R_radar_;



  // create matrix for cross correlation Tc
  Eigen::MatrixXd Tc = Eigen::MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  Eigen::MatrixXd K = Tc * S.inverse();

  // residual
  Eigen::VectorXd ZRealdiff = z - z_pred;

  // angle normalization
  while (ZRealdiff(1) > M_PI) ZRealdiff(1) -= 2. * M_PI;
  while (ZRealdiff(1) < -M_PI) ZRealdiff(1) += 2. * M_PI;

  // update state mean and covariance matrix
  x_ += K * ZRealdiff;
  P_ -= K * S * K.transpose();

  // Calculate Lidar NIS
  NIS_radar_ = ZRealdiff.transpose() * S.inverse() * ZRealdiff;
}


void UKF::PredictSigmaPoints(double delta_t) {
  // create augmented mean vector
  Eigen::VectorXd x_aug = Eigen::VectorXd(n_aug_);

  // create augmented state covariance
  Eigen::MatrixXd P_aug = Eigen::MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);

  // create augmented sigma point matrix
  Eigen::MatrixXd Xsig_aug = Eigen::MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);

  // create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  // create square root matrix
  Eigen::MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int t = 0; t < n_aug_; t++) {
    Xsig_aug.col(t + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(t);
    Xsig_aug.col(t + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(t);
  }

  // predict sigma points
  for (int i = 0; i < (2 * n_aug_ + 1); ++i) {
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p {};
    double py_p {};

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + (v/yawd) * (sin(yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + (v/yawd) * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;


    // add noise
    px_p = px_p + 0.5 * nu_a * powf(delta_t, 2) * cos(yaw);
    py_p = py_p + 0.5 * nu_a * powf(delta_t, 2) * sin(yaw);
    v_p = v_p + nu_a * delta_t;
    yaw_p = yaw_p + 0.5 * nu_yawdd * powf(delta_t, 2);
    yawd_p = yawd_p + nu_yawdd * delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

void UKF::PredictStateMeanAndCovariance() {
  // predicted state mean
  x_.noalias() = Xsig_pred_ * weights_;

  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}
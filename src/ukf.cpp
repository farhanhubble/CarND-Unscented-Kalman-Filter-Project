#include "ukf.h"
#include "Eigen/Dense"
#include "tools.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  // State initialization flag.
  is_initialized_ = false;

  // Time stamp of last sensor package received.
  long long time_us_ = 0;

  // Number of state variables.
  n_x_ = x_.rows();

  // Number of states in the augmented state vector.
  n_aug_ = n_x_ + 2;

  // Number of sigma points.
  n_sigma_ = 2 * n_aug_ + 1;

  // Spread parameter for sigma points generation.
  lambda_ = 3 - n_x_;

  // Sigma point weights.
  weights_ = VectorXd(n_sigma_);
  weights_.head(1) << lambda_ / (lambda_ + n_aug_);
  weights_.tail(n_sigma_-1) << 1 / (2*(lambda_ + n_aug_));

  // Predicted sigma point matrix.
  Xsig_pred_ = MatrixXd(n_x_,n_sigma_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR && 
    use_radar_ == false) {
      return;
    }
  if(meas_package.sensor_type_ == MeasurementPackage::LASER && 
    use_laser_ == false) {
      return;
    }


  if(is_initialized_ == false){ 
    initialize_state(meas_package);
  }
  else{
    double delta_t_secs = (meas_package.timestamp_ - time_us_) * 1e-6;
    Prediction(delta_t_secs);

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){ 
      UpdateRadar(meas_package); 
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      UpdateLidar(meas_package);
    }
  }


  time_us_ = meas_package.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // Construct augmented mean state vector by appending mean values of tangential.
  // and normal acceleration (=0).
  Eigen::VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  for(int i=n_x_; i<n_aug_; i++){
    x_aug(i) = 0;
  }

  // Construct the process noise covariance matrix;
  Eigen::MatrixXd Q;
  Q << std_a_*std_a_, 0,
       0, std_yawdd_*std_yawdd_;
  // Construct augmented covariance matrix P_aug.
  Eigen::MatrixXd P_aug = MatrixXd::Zero(n_aug_,n_aug_);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug.bottomRightCorner(n_aug_-n_x_,n_aug_-n_x_) = Q;

  // Construct the square root matrix of P_aug.
  Eigen::MatrixXd A = P_aug.llt().matrixL();
  A = A * sqrt(lambda_+n_aug_);

  // Construct the sigma point matrix X_sigma.
  MatrixXd X_sigma = MatrixXd(n_aug_,n_aug_);
  // Fill first column with the augmented mean state vector x_aug.
  X_sigma.col(0) =  x_aug;
  // Fill second to n_aug columns  with x_aug+A;
  X_sigma.block<7,7>(0,1) = A.colwise() + x_aug;
  // Fill n_aug+1 to 2*n_aug+1 columns with x_aug-A;
  X_sigma.block<7,7>(0,8) = (-1*A).colwise() + x_aug;

  /* ACTUAL PREDICTION STEP. */

  // Generate output sigma points by passing each column of the generated sigma
  // points matrix through the (non-linear) process function f(.)
  for(int i=0; i<n_sigma_; i++){
    double px      = X_sigma(0,i); // Cartesian x coordinate.
    double py      = X_sigma(1,i); // Cartesian y coordinate.
    double v       = X_sigma(2,i); // Tangential velocity.
    double psi     = X_sigma(3,i); // Turn angle.
    double psi_dot = X_sigma(4,i); // Turn rate (angular velocity).
    double nu_tang = X_sigma(5,i); // Tangential acceleration.
    double nu_norm = X_sigma(6,i); // Angular acceleration.

    // Compute noise vector.
    MatrixXd noise(n_x_,1);
    noise << // Effect of tangential acceleration on x position.
            0.5 * nu_tang * cos(psi) * delta_t * delta_t, 
            // Effect of tangential acceleration on y position.
            0.5 * nu_tang * sin(psi) * delta_t * delta_t,
            // Effect of tangential acceleration tangential velocity.
            nu_tang * delta_t,
            // Effect of angular acceleration on turn angle.
            0.5 * nu_norm * delta_t * delta_t;
            // Effect of angular acceleration on angular velocity.
            nu_norm * delta_t;

    // Compute update vector.
                            
    MatrixXd update(n_x_,1);
    if(psi_dot == 0){
        update << v*cos(psi)*delta_t,
               v*sin(psi)*delta_t,
               0,
               psi_dot*delta_t,
               0;
    }
    else{
        update <<  v/psi_dot*(sin(psi+psi_dot*delta_t) - sin(psi)),
               v/psi_dot*(-cos(psi+psi_dot*delta_t) + cos(psi)),
               0,
               psi_dot*delta_t,
               0;
    }

    // Fill in the Sigma point prediction matrix.
    Xsig_pred_.col(i) = X_sigma.col(i) + update + noise;

    // Calculate mean state vector as weighted mean of columns of the 
    // **predicted** sigma points matrix.
    x_ = ( Xsig_pred_.array().rowwise() * weights_.array().transpose() ).rowwise().sum();
    // Normalize angles.
    x_.row(3).unaryExpr(&Tools::NormalizeAngle);


    // Calculate state covariance matrix as weighted covariance of the rows
    // of the **predicted** sigma points matrix.
    MatrixXd Error = Xsig_pred_.array().colwise() - x_.array();
    MatrixXd WeightedError = Error.array().rowwise() * weights_.array().transpose();
    P_ = WeightedError * Error.transpose();
  }
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // Reuse the **predicted** sigma points from the prediction step to get the 
  // measurement sigma points Z_sigma. Radar has 3-dimensional measurement and
  // measurement noise is additive so no augmentation is performed and the shape
  // of the sigma point matrix will be (3,2*3+1)

  // Actual measurement vector z.
  int n_z_ = meas_package.raw_measurements_.size();
  VectorXd z = VectorXd(n_z_);
  
  MatrixXd Zsig_pred = MatrixXd::Zero(n_z_,n_sigma_);

  // Fill in the columns of the measurement sigma points matrix.
  for(int i=0; i<n_sigma_; i++){
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double psi = Xsig_pred_(3,i);
    double psi_dot = Xsig_pred_(4,i);

    // Calculate prediction space variables.
    
    double rho = sqrt(px*px+py*py);
    // Handle division by zero
    if(rho == 0){
      rho = 1e-2;
    }

    double phi = atan2(py,px);
    assert((phi <= M_PI) && (phi >= -M_PI)); 

    double rho_dot = (v*cos(psi)*px + v*sin(psi)*py) / rho;

    Zsig_pred.col(i) << rho, phi, rho_dot; 
  }

  // Calculate mean predicted measument vector z as weighted mean of the
  // prediction space sigma point matrix Zsig_pred. 
  VectorXd z_pred = VectorXd(n_z_);
  z_pred = (Zsig_pred.array().rowwise() * weights_.array().transpose()).rowwise().sum();

  //Calculate measurement covariance matrix S as weighted covariance of the 
  // prediction space sigma point matrix Zsig_pred.

  // Compute the measurement noise covariancce matrix R.
  MatrixXd R = MatrixXd::Zero(n_z_,n_z_);
  R.diagonal() << std_radr_*std_radr_, std_radphi_*std_radphi_,  std_radrd_*std_radrd_;

  // Compute the covariance matrix of predicted measurement.
  MatrixXd S = MatrixXd(n_z_,n_z_);
  MatrixXd Error = Zsig_pred.array().colwise() - z_pred.array();
  Error.row(1).unaryExpr(&Tools::NormalizeAngle);
  MatrixXd WeightedError = Error.array().rowwise() * weights_.array().transpose();
  S = (WeightedError * Error.transpose()) + R;

  /*------------------------ APPLY ACTUAL UPDATE-------------------------------- */
  
  // Calculate the cross correlation matrix Tc of Xsig_pred and Zsig_pred;
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  MatrixXd Error_X = Xsig_pred_.array().colwise() - x_.array();
  Error_X.row(3).unaryExpr(&Tools::NormalizeAngle);
  MatrixXd WeightedError_X = Error_X.array().rowwise() * weights_.array().transpose();
  
  MatrixXd Error_Z = Zsig_pred.array().colwise() - z_pred.array();
  Error_Z.row(1).unaryExpr(&Tools::NormalizeAngle);

  Tc = WeightedError_X * Error_Z.transpose();


  // Find kalman gain matrix K.
  MatrixXd K = Tc * S.inverse(); 

  //update state mean and covariance matrix
  VectorXd y = z-z_pred;
  y.row(1).unaryExpr(&Tools::NormalizeAngle);
  x_ = x_ + K * y;
  P_ = P_ - K * S * K.transpose();
}

void UKF::initialize_state(MeasurementPackage meas_package){
  x_ << 0,0,0,0,0;
  P_ << 1,0,0,0,0,
        0,1,0,0,0,
        0,0,1,0,0,
        0,0,0,1,0,
        0,0,0,0,1;

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    double rho = meas_package.raw_measurements_[0];
    double phi = meas_package.raw_measurements_[1];
    x_(0) = rho*cos(phi);
    x_(1) = rho*sin(phi);

  }

  else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    x_(0) = meas_package.raw_measurements_[0];
    x_(1) = meas_package.raw_measurements_[1];
  }
}
#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // Constructor so not initialised yet
  HasInitialMeasurement_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  UseLaser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  UseRadar_ = true;

  // state vector length
  NumX_ = 5;

  // augmentation vector length
  NumAug_ = 7;

  // lambda value
  Lambda_ = 3 - NumAug_;

  // initial state vector
  X_ = VectorXd(NumX_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(NumX_, NumX_);

  // initial weights
  Weights_ = VectorXd((2 * NumAug_) + 1);

  // initial sigma predicition matrix
  XSigPred_ = MatrixXd(NumX_, 2 * NumAug_ + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // Value from "Augmentation Assignment 2" lecture
  StdA_ = 0.75;
  //std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // Value from "Augmentation Assignment 2" lecture
  StdYawdd_ = 0.75;

  // Laser measurement noise standard deviation position1 in m
  StdLasPx_ = 0.1;

  // Laser measurement noise standard deviation position2 in m
  StdLasPy_ = 0.1;

  // Radar measurement noise standard deviation radius in m
  StdRadR_ = 0.2;

  // Radar measurement noise standard deviation angle in rad
  StdRadPhi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  StdRadRd_ = 0.75;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  
  // Check the measurement type is from a valid source.
  assert((meas_package.sensor_type_ == MeasurementPackage::RADAR) ||
         (meas_package.sensor_type_ == MeasurementPackage::LASER));

  // If the state vector has not been initialised - use the first measurement
  if (!HasInitialMeasurement_)
  {
    VectorXd InitialState(NumX_);

    // first measurement
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      // Convert radar from polar to cartesian coordinates and initialize state.

      // Used for positional calculation
      float Rho     = meas_package.raw_measurements_[0];
      // Used for velocity calculation
      float RhoDot  = meas_package.raw_measurements_[2];
      // Angle for polar coordinates
      float Theta   = meas_package.raw_measurements_[1];

      float PositionX = Rho * cos(Theta);
      float PositionY = Rho * sin(Theta);
      float VelocityX = RhoDot * cos(Theta);
      float VelocityY = RhoDot * sin(Theta);

      InitialState << PositionX, PositionY, VelocityX, VelocityY, 0;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      // Already in cartesian coordinates - so don't convert
      float PositionX  = meas_package.raw_measurements_[0];
      float PositionY  = meas_package.raw_measurements_[1];

      InitialState << PositionX, PositionY, 0, 0, 0;
    }

    // Initialise the kalman filter state
    Feed(InitialState, meas_package.timestamp_);
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  Prediction(meas_package.timestamp_);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates
    UpdateRadar(meas_package);
  } 
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    // Laser updates
    UpdateLidar(meas_package);
  }
}

/**
 * @brief Feed the kalman filter with the first measurement
 *
 * @param initialState The first measurement
 * @param timestamp The time of the first measurement
 */
void UKF::Feed(const Eigen::VectorXd &InitialState,
               long long Timestamp)
{
  // Set the state vector
  X_ = InitialState;

  HasInitialMeasurement_ = true;
  PreviousTimestamp_     = Timestamp;
}

/**
 * Create the mean augmented state vector
 */
void UKF::CreateAugmentedMeanState(VectorXd &XAug)
{
  // Clear the augmented state
  XAug.fill(0);

  // Update te augmented state
  XAug.head(5) = X_;
}

/**
 * Create the augmented covariance matrix
 */
void UKF::CreateAugmentedCovarience(MatrixXd &PAug)
{
  // Clear the covariances matrix
  PAug.fill(0.0);

  // Update the augmented covariance matrix
  PAug.topLeftCorner(5,5) = P_;
  PAug(5,5)               = StdA_ * StdA_;
  PAug(6,6)               = StdYawdd_ * StdYawdd_;
}

/**
 * Create sqaure root matrix
 */
void UKF::CreateSqrtMatrix(MatrixXd &Sqrt,
                           MatrixXd &PAug)
{
  Sqrt = PAug.llt().matrixL();
}

/**
 * Create augmented sigma points
 */
void UKF::CreateAugmentedSigmaPoints(MatrixXd &XSigAug,
                                     VectorXd &XAug,
                                     MatrixXd &Sqrt)
{
  int i = 0;

  // Clear augmented sigma matrix
  XSigAug.fill(0.0);

  // Determine augmented sigma points
  XSigAug.col(0)  = XAug;
  for (i = 0; i < NumAug_; i++)
  {
    XSigAug.col(i + 1)           = XAug + sqrt(Lambda_ + NumAug_) * Sqrt.col(i);
    XSigAug.col(i + 1 + NumAug_) = XAug - sqrt(Lambda_ + NumAug_) * Sqrt.col(i);
  }
}

/**
 * Create predicited sigma points
 */
void UKF::CreatePredictedSigmaPoints(MatrixXd &XSigAug,
                                     double TimeStep)
{
  int i = 0;

  int SigmaPoints = XSigAug.size();

  // Predict sigma points
  for (i = 0; i < (2 * NumAug_) + 1; i++)
  {
    // Extract values for better readability
    double PositionX = XSigAug(0,i);
    double PositionY = XSigAug(1,i);
    double Velocity  = XSigAug(2,i);
    double Yaw       = XSigAug(3,i);
    double Yawd      = XSigAug(4,i);
    double NuYaw     = XSigAug(5,i);
    double NuYawdd   = XSigAug(6,i);

    // Predicited values
    double PredictPositionX, PredictPositionY;
    double PredictVelocity;
    double PredictYaw, PredictYawd;

    // Predicted Position Values
    // Note: Calculation is different if Yawd is 0 (possible divide by zero)
    if (fabs(Yawd) > numeric_limits<float>::epsilon())
    {
      PredictPositionX = PositionX + Velocity/Yawd * (sin(Yaw + (Yawd * TimeStep)) - sin(Yaw));
      PredictPositionY = PositionY + Velocity/Yawd * (cos(Yaw) - cos(Yaw + (Yawd * TimeStep)));
    }
    else
    {
      PredictPositionX = PositionX + (Velocity * TimeStep * cos(Yaw));
      PredictPositionY = PositionY + (Velocity * TimeStep * sin(Yaw));
    }

    // Predict Velocity
    PredictVelocity = Velocity;

    // Predict Yaw
    PredictYaw      = Yaw + (Yawd * TimeStep);
    PredictYawd     = Yawd;

    // Add noise to predictions
    PredictPositionX = PredictPositionX + (0.5 * NuYaw * TimeStep * TimeStep * cos(Yaw));
    PredictPositionY = PredictPositionY + (0.5 * NuYaw * TimeStep * TimeStep * sin(Yaw));
    PredictVelocity  = PredictVelocity  + (NuYaw * TimeStep);
    PredictYaw       = PredictYaw  + (0.5 * NuYawdd * TimeStep * TimeStep);
    PredictYawd      = PredictYawd + (NuYawdd * TimeStep);

    // Write predicted sigma point into right column
    XSigPred_(0,i) = PredictPositionX;
    XSigPred_(1,i) = PredictPositionY;
    XSigPred_(2,i) = PredictVelocity;
    XSigPred_(3,i) = PredictYaw;
    XSigPred_(4,i) = PredictYawd;
  }
}

/**
 * Set Weights
 */
void UKF::SetWeights()
{
  int i = 0;

  // set weights
  Weights_(0) = Lambda_ / (Lambda_ + NumAug_);

  // (2n + 1) weights
  for (i = 1; i < (2 * NumAug_) + 1; i++)
  {  
    Weights_(i) = 0.5 / (NumAug_ + Lambda_);
  }
}

/**
 * Predict State Mean
 */
void UKF::PredictStateMean()
{
  int i = 0;

  // Clear state vector
  X_.fill(0.0);

  // Iterate over sigma points
  for (i = 0; i < (2 * NumAug_) + 1; i++)
  {  
    // Update state vector
    X_ = X_ + Weights_(i) * XSigPred_.col(i);
  }
}

/**
 * Predict State Covariance
 */
void UKF::PredictStateCovariance()
{
  int i = 0;

  // Clear covariance matrix
  P_.fill(0.0);

  // Iterate over sigma points
  for (i = 0; i < (2 * NumAug_) + 1; i++)
  {  
    // State difference
    VectorXd XDiff = XSigPred_.col(i) - X_;

    // Angle normalization
    XDiff(3) = tools.NormalizeAngle(XDiff(3));

    // Update covariance matrix
    P_ = P_ + (Weights_(i) * XDiff * XDiff.transpose());
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param timestamp The timestamp to predict to
 */
void UKF::Prediction(long long timestamp)
{
  //create augmented mean vector
  VectorXd XAug = VectorXd(NumAug_);

  //create augmented state covariance
  MatrixXd PAug = MatrixXd(NumAug_, NumAug_);

  //create sigma point matrix
  MatrixXd XSigAug = MatrixXd(NumAug_, (2 * NumAug_) + 1);

  // Square root matrix
  MatrixXd Sqrt;

  // Calculate the time offset for the prediction
  double PredictTimeStep =  ((float)(timestamp - PreviousTimestamp_))/1e6;
  PreviousTimestamp_ = timestamp;
 
  // Create augmented mean state
  CreateAugmentedMeanState(XAug);

  // Create augmented covariance matrix
  CreateAugmentedCovarience(PAug);

  // Populate square root matrix
  CreateSqrtMatrix(Sqrt, PAug);

  // Create augmented sigma points
  CreateAugmentedSigmaPoints(XSigAug, XAug, Sqrt);

  // Predict Sigma Points
  CreatePredictedSigmaPoints(XSigAug, PredictTimeStep);

  /****************************************************
   *                                                  *
   *           Predict Mean and covariance            *
   *                                                  *
   ****************************************************/

  // Set Weights
  SetWeights();

  // Predicted state mean
  PredictStateMean();

  // Predicted state covariance matrix
  PredictStateCovariance();
}

/**
 * Update the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  int i = 0;

  // Set measurement dimension, lidar can measure px, py
  int NumMeasurementDims = 2;

  // Create matrix for sigma points in measurement space
  MatrixXd ZSig = MatrixXd(NumMeasurementDims, (2 * NumAug_) + 1);

  // Prediction state
  VectorXd ZPred = VectorXd(NumMeasurementDims);

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(NumMeasurementDims, NumMeasurementDims);

  // Cross correlation matrix Tc
  MatrixXd Tc = MatrixXd(NumX_, NumMeasurementDims);

  // Noise covariance matrix
  MatrixXd R = MatrixXd(NumMeasurementDims, NumMeasurementDims);

  // Kalman gain
  MatrixXd K; 

  // Residual
  VectorXd Residual;

  // R values
  double StdLasPx2, StdLasPy2;

  // NIS
  double NISLidar;

  // Clear all matrices
  S.fill(0.0);
  Tc.fill(0.0);
  ZSig.fill(0.0);
  ZPred.fill(0.0);

  // Transform sigma points into measurement space
  // (2n+1) simga points
  for (i = 0; i < (2 * NumAug_) + 1; i++) 
  {  
    // Extract values for better readibility
    double PositionX = XSigPred_(0,i);
    double PositionY = XSigPred_(1,i);

    // Measurement model
    ZSig(0,i) = PositionX;
    ZSig(1,i) = PositionY;
  }

  // Mean predicted measurement
  for (i = 0; i < (2 * NumAug_) + 1; i++)
  {
    ZPred = ZPred + (Weights_(i) * ZSig.col(i));
  }
  
  // 2n + 1 simga points
  for (i = 0; i < (2 * NumAug_) + 1; i++) 
  {  
    // Residual
    VectorXd _Residual = ZSig.col(i) - ZPred;

    // State difference
    VectorXd XDiff = XSigPred_.col(i) - X_;

    // Angle normalization
    _Residual(1) = tools.NormalizeAngle(_Residual(1));

    // Angle normalization
    XDiff(3) = tools.NormalizeAngle(XDiff(3));

    Tc = Tc + Weights_(i) * XDiff * _Residual.transpose();
    S  = S  + Weights_(i) * _Residual * _Residual.transpose();
  }

  // Add measurement noise covariance matrix
  StdLasPx2 = (StdLasPx_ * StdLasPx_);
  StdLasPy2 = (StdLasPy_ * StdLasPy_);

  R <<  StdLasPx2, 0,
        0,         StdLasPy2;
  
  // Covariance matrix
  S = S + R;

  /****************************************************
   *                                                  *
   *                   Update State                   *
   *                                                  *
   ****************************************************/

  // Kalman gain K;
  K = Tc * S.inverse();

  // Residual
  Residual = meas_package.raw_measurements_ - ZPred;

  // Angle normalization
  Residual(1) = tools.NormalizeAngle(Residual(1));

  // Update state mean and covariance matrix
  X_ = X_ + (K * Residual);
  P_ = P_ - (K * S * K.transpose());

  // Calculate NIS
  NISLidar = Residual.transpose() * S.inverse() * Residual;
  cout << "NIS_LIDAR " << NISLidar << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  int i = 0;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int NumMeasurementDims = 3;
 
  // Create matrix for sigma points in measurement space
  MatrixXd ZSig = MatrixXd(NumMeasurementDims, (2 * NumAug_) + 1);

  // Prediction state
  VectorXd ZPred = VectorXd(NumMeasurementDims);

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(NumMeasurementDims, NumMeasurementDims);

  // Cross correlation matrix Tc
  MatrixXd Tc = MatrixXd(NumX_, NumMeasurementDims);

  // Noise covariance matrix
  MatrixXd R = MatrixXd(NumMeasurementDims, NumMeasurementDims);

  // Kalman gain
  MatrixXd K; 

  // Residual
  VectorXd Residual;

  // R values
  double StdRadR2, StdRadPhi2, StdRadRd2;

  // NIS
  double NISRadar;

  // Clear all matrices
  S.fill(0.0);
  Tc.fill(0.0);
  ZSig.fill(0.0);
  ZPred.fill(0.0);

  //transform sigma points into measurement space
  //2n+1 simga points
  for (i = 0; i < (2 * NumAug_) + 1; i++)
  {  
    // extract values for better readibility
    double PositionX = XSigPred_(0,i);
    double PositionY = XSigPred_(1,i);
    double Velocity  = XSigPred_(2,i);
    double Yaw       = XSigPred_(3,i);

    double Velocity1 = cos(Yaw) * Velocity;
    double Velocity2 = sin(Yaw) * Velocity;

    // measurement model
    // rho
    ZSig(0,i) = sqrt((PositionX * PositionX) + (PositionY * PositionY));
    // phi
    ZSig(1,i) = atan2(PositionY, PositionX);
    // rho_dot
    ZSig(2,i) = ((PositionX * Velocity1) + (PositionY * Velocity2)) / ZSig(0,i);
  }

  for (i = 0; i < (2 * NumAug_) + 1; i++)
  {
    ZPred = ZPred + Weights_(i) * ZSig.col(i);
  }

  for (int i = 0; i < (2 * NumAug_) + 1; i++) 
  {  
    //2n+1 simga points
    //residual
    VectorXd _Residual = ZSig.col(i) - ZPred;

    //angle normalization
    _Residual(1) = tools.NormalizeAngle(_Residual(1));

    // state difference
    VectorXd XDiff = XSigPred_.col(i) - X_;
    
    //angle normalization
    XDiff(3) = tools.NormalizeAngle(XDiff(3));

    Tc = Tc + Weights_(i) * XDiff * _Residual.transpose();
    S  = S + Weights_(i) * _Residual * _Residual.transpose();
  }

  //add measurement noise covariance matrix
  StdRadR2   = (StdRadR_ * StdRadR_);
  StdRadPhi2 = (StdRadPhi_ * StdRadPhi_);
  StdRadRd2  = (StdRadRd_ * StdRadRd_);
  R <<  StdRadR2, 0,          0,
        0,        StdRadPhi2, 0,
        0,        0,          StdRadRd2;


  // Covariance matrix
  S = S + R;

  /****************************************************
   *                                                  *
   *                   Update State                   *
   *                                                  *
   ****************************************************/

  //Kalman gain K;
  K = Tc * S.inverse();

  //residual
  Residual = meas_package.raw_measurements_ - ZPred;

  //angle normalization
  Residual(1) = tools.NormalizeAngle(Residual(1));

  //update state mean and covariance matrix
  X_ = X_ + (K * Residual);
  P_ = P_ - (K * S * K.transpose());

  // Calculate NIS
  NISRadar = Residual.transpose() * S.inverse() * Residual;
  cout << "NIS_RADAR " << NISRadar << endl;
}

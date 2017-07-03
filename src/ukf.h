#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include "tools.h"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool HasInitialMeasurement_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool UseLaser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool UseRadar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd X_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd XSigPred_;

  ///* time when the state is true, in us
  long long PreviousTimestamp_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double StdA_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double StdYawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double StdLasPx_;

  ///* Laser measurement noise standard deviation position2 in m
  double StdLasPy_;

  ///* Radar measurement noise standard deviation radius in m
  double StdRadR_;

  ///* Radar measurement noise standard deviation angle in rad
  double StdRadPhi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double StdRadRd_ ;

  ///* Weights of sigma points
  VectorXd Weights_;

  ///* State dimension
  int NumX_;

  ///* Augmented state dimension
  int NumAug_;

  ///* Sigma point spreading parameter
  double Lambda_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * @brief Feed the kalman filter with the first measurement
   *
   * @param initialState The first measurement
   * @param timestamp The time of the first measurement
   */
  void Feed(const Eigen::VectorXd &initialState,
            long long timestamp);

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);




  void CreateAugmentedMeanState(VectorXd &XAug);

  void CreateAugmentedCovarience(MatrixXd &PAug);

  void CreateSqrtMatrix(MatrixXd &Sqrt,
                        MatrixXd &PAug);

  void CreateAugmentedSigmaPoints(MatrixXd &XSigAug,
                                  VectorXd &XAug,
                                  MatrixXd &Sqrt);

  void CreatePredictedSigmaPoints(MatrixXd &XSigAug,
                                  double TimeStep);

  void SetWeights();

  void PredictStateMean();

  void PredictStateCovariance();

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(long long timestamp);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

private:
  Tools tools;
};

#endif /* UKF_H */

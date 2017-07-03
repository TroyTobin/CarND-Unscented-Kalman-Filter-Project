#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
  // Create the RMSE output vector and initialize
  uint32_t i;
  VectorXd rmse(4);

  rmse << 0,0,0,0;

  // Sanity check the inputs that we have the same number of 
  // estimations as we do reference data
  if (estimations.size() != ground_truth.size())
  {
    cerr << "CalculateRMSE: estimate and ground_truth size mismatch ";
    cerr << "(" << estimations.size() << ", " << ground_truth.size() << ")" << endl;
    goto Error;
  }

  // if there are no estimations nothing to do
  if (estimations.size() == 0)
  {
    cerr << "CalculateRMSE: No estimations ";
    cerr << "(" << estimations.size() << ")" << endl;
    goto Error;
  }

  // Sum the squared residuals for each sample
  for(i = 0; i < estimations.size(); i++)
  {
    // Sanity check the the samples are the correct sizes
    VectorXd est_sample = estimations[i];
    VectorXd gt_sample  = ground_truth[i];

    if ((est_sample.size() != 4) ||
      (gt_sample.size()  != 4))
    {
      cerr << "CalculateRMSE: Invalid estimate and ground_truth sample size ";
      cerr << "(" << est_sample.size() << ", " << gt_sample.size() << ") != " << 4 << endl;
      goto Error;     
    }

    VectorXd residual = est_sample - gt_sample;

    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

Error:
  //return the result
  return rmse;
}



/**
 * A helper method to normalize an angle to [-pi, pi]
 */
double Tools::NormalizeAngle(double angle)
{
  // Angle normalization
  while (angle > M_PI)
  {
    angle -= 2.*M_PI;
  }
  while (angle < -M_PI)
  {
    angle += 2.*M_PI;
  }
  
  return angle;
}
#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
}

double Tools::NormalizeAngle(const double theta){
  double result = theta;
  while(result < -M_PI) result += 2*M_PI;
  while(result >  M_PI) result -= 2*M_PI;

  return result;
}
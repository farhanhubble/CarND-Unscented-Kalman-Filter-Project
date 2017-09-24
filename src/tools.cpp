#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth){
  assert( (estimations.size() !=0 ) && (estimations.size() == ground_truth.size()) );
   
  VectorXd RMSE = VectorXd::Zero(estimations[0].size());
  for(int i=0; i<estimations.size(); i++){
    VectorXd error = estimations[i]-ground_truth[i];
    error = error.array() * error.array();
    RMSE += error;
  }

  RMSE = RMSE / estimations.size();
  RMSE = RMSE.array().sqrt();
  return RMSE;
}

double Tools::NormalizeAngle(const double theta){
  double result = theta;
  while(result < -M_PI) result += 2*M_PI;
  while(result >  M_PI) result -= 2*M_PI;

  return result;
}
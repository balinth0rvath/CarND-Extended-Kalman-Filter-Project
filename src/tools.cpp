#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// Must be similar size, must contain elements
	if (!estimations.size() || estimations.size()!=ground_truth.size())
		return rmse;

	VectorXd sum(4);
	sum << 0,0,0,0;

	// TODO: optimize code eliminating temporary varibles 
	for(unsigned int i=0; i < estimations.size(); ++i)
	{
		VectorXd sb = estimations[i] - ground_truth[i];
		sum += (VectorXd)(sb.array() * sb.array());
	}	
	rmse = sum/estimations.size();
	rmse = rmse.array().sqrt();

	return rmse;	
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
	MatrixXd Hj(3,4);
	
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	if (px == 0 && py == 0)
		return Hj;

	double p_pow = pow(px,2) + pow(py,2);
	double pv_mult_1 = py * (vx * py - vy * px);
	double pv_mult_2 = px * (vy * px - vx * py);

	Hj <<	px/sqrt(p_pow),						py/sqrt(p_pow),						0,							0,
				-py/p_pow,								px/p_pow,									0,							0,	
				pv_mult_1/pow(p_pow,1.5), pv_mult_2/pow(p_pow,1.5), px/sqrt(p_pow), py/sqrt(p_pow);

	return Hj;
}


Eigen::VectorXd Tools::ConvertToPolar(const Eigen::VectorXd& cartesian)
{
	VectorXd polar = VectorXd(3);

	double rho, phi, rho_dot;
	double p_x = cartesian[0];
	double p_y = cartesian[1];
	double v_x = cartesian[2];
	double v_y = cartesian[3];

	rho = sqrt(p_x * p_x+ p_y * p_y);
	phi = atan2(p_y, p_x);

	rho = (rho < 0.0001) ? 0.0001 : rho;
	rho_dot = (p_x * v_x + p_y * v_y) / rho;
	polar << rho, phi, rho_dot;
				
	return polar;
}

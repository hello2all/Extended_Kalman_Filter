#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::atan;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_laser_in, MatrixXd &R_radar_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_laser_ = R_laser_in;
  R_radar_ = R_radar_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  	MatrixXd Hj = tools.CalculateJacobian(x_);
	VectorXd z_pred = this->Cartesian2Polar(x_);
	VectorXd y = z - z_pred;
	MatrixXd Ht = Hj.transpose();
	MatrixXd S = Hj * P_ * Ht + R_radar_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj) * P_;
}

VectorXd KalmanFilter::Cartesian2Polar(const VectorXd &x){
	float px = x(0);
	float py = x(1);
	float vx = x(2);
	float vy = x(3);

	float rho = sqrt(px * px + py * py);
	float phi = atan(px / py);
	float rho_dot = (px * vx + py * vy) / sqrt(px * px + py * py);

	Eigen::VectorXd h_x_ = VectorXd(3);
	h_x_ << rho, phi, rho_dot;

	return h_x_;
}
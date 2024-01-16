#pragma once
#include<Eigen/Dense>
#include<Eigen/Core>
#include<vector>

class solution
{
public:
	void qp_lagrange(Eigen::MatrixXd& H, Eigen::MatrixXd& c, Eigen::MatrixXd& A, Eigen::MatrixXd& b,
		Eigen::MatrixXd& x, Eigen::MatrixXd& lambda, const int& dim, const int& m);

	std::vector<Eigen::MatrixXd> eigen_elimination(Eigen::MatrixXd A, Eigen::MatrixXd b) {
		Eigen::HouseholderQR<Eigen::MatrixXd> qr;
		qr.compute(A);
		Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
		Eigen::MatrixXd Q = qr.householderQ();

		//Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(R);
		//int rank = lu_decomp.rank();

		//if (!rank == A.rows()) {
		//	printf("##########\n");
		//}

		int cols = R.cols();
		Eigen::MatrixXd retb = Q.inverse()*b;
		Eigen::MatrixXd retbb = retb.block(0, 0, 28, 1);

		return { R.block(0,0,28,cols),retbb };
	}

};


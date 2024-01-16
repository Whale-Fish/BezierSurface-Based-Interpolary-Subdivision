#pragma once

#include "eigenSolution.h"
#include<iostream>
#include <fstream>
extern bool testFlg;

void solution::qp_lagrange(Eigen::MatrixXd& H, Eigen::MatrixXd& c, Eigen::MatrixXd& A, Eigen::MatrixXd& b,
	Eigen::MatrixXd& x, Eigen::MatrixXd& lambda, const int& dim, const int& m) {

	Eigen::MatrixXd G(dim, dim);
	Eigen::MatrixXd B(m, dim);
	Eigen::MatrixXd C(m, m);

	Eigen::MatrixXd invH = H.inverse();
	Eigen::MatrixXd transA = A.transpose();

	Eigen::MatrixXd tmp = ((A * invH * transA));

	double d = tmp.determinant();

	Eigen::MatrixXd invTmp = tmp.inverse();

	G = invH - invH * transA * invTmp * A * invH;
	B = (A * invH * transA).inverse() * A * invH;
	//C = -1 * (A * invH * transA).inverse();

	if (testFlg) {
		//std::cout << "===========================" << std::endl;
		//std::cout << invH << std::endl;
		//std::cout << "===========================" << std::endl;
		//std::cout << "===========================" << std::endl;
		//std::cout << tmp << std::endl;
		//std::cout << "===========================" << std::endl;
		//std::cout << "===========================" << std::endl;
		//std::cout << transA << std::endl;
		//std::cout << "===========================" << std::endl;
		//std::cout << "===========================" << std::endl;
		//std::cout << G << std::endl;
		//std::cout << "===========================" << std::endl;

		//std::ofstream outfile("matrix.txt"); // ��������ļ���
		//outfile << tmp.format(Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "\n", "", "", "", "")); // ��������txt��ʽд���ļ�
		//outfile.close(); // �ر��ļ���
	}

	x = B.transpose() * b - G * c;
	//lambda = B * c - C * b;

}

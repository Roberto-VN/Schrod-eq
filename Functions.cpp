#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <armadillo>
#include <tuple>
#include <complex>

#include "Functions.hpp"

// Define Create_Matrix function

std::tuple<arma::sp_cx_mat, arma::sp_cx_mat> Create_Matrix(double h, double dt, arma::mat V)
{
// Define M - 2 & (M - 2)^2
	
	int M = 1. / h;
	int M_1 = M - 2; // Dimension of the submatrices
	int M_2 = M_1 * M_1; // Dimension of the matrices

// Define r

	arma::cx_double r(0.0, dt / (2 * h * h));
	
// Define the submatrices

	arma::cx_mat submat1(M_1, M_1, arma::fill::zeros);
	arma::cx_mat submat2(M_1, M_1, arma::fill::zeros);

	// Define the diagonal & subdiagonals of the submatrices
	
	arma::cx_vec diag(M_1, arma::fill::zeros);
	diag.fill(arma::cx_double(1.0, 0.0));
	diag *= r;

	arma::cx_vec subdiag(M_1 - 1, arma::fill::zeros);
	subdiag.fill(arma::cx_double(1.0, 0.0));
	subdiag *= r;

	// Create the submatrices

	submat1 += arma::diagmat(subdiag, -1);
	submat1 += arma::diagmat(subdiag, +1);
	submat2 += arma::diagmat(diag);

// Define the matrices A & B

	// Create a matrix that that translates the indices (i, j) into k

	arma::imat K(M_1, M_1, arma::fill::zeros);

	int k_index = 0;

	for (int i = 0; i < M_1; i++)
	{
		for (int j = 0; j < M_1; j++)
		{
			K(i, j) = k_index;
			k_index += 1;
		}
	}
	
	// Create the Matrices 

	arma::sp_cx_mat A(M_2, M_2);

	for (int i = 0; i < M_2; i += M_1) 
	{
		A.submat(i, i, i + M_1 - 1, i + M_1 - 1) = -1. * submat1;
	}

	for (int i = 0; i < M_2 - M_1; i += M_1)
	{
		A.submat(i, i + M_1, i + M_1 - 1, i + 2 * M_1 - 1) = -1. * submat2;
		A.submat(i + M_1, i, i + 2 * M_1 - 1, i + M_1 - 1) = -1. * submat2;
	}

	arma::sp_cx_mat B = -1. * A;

	// Define the diagonals of the matrices

	arma::cx_vec a(M_2);
	arma::cx_vec b(M_2);

	for (int i = 0; i < M_1; i++)
	{
		for (int j = 0; j < M_1; j++)
		{
			a(K(i, j)) = arma::cx_double(1.0, dt * V(i, j) / 2) + 4. * r;
			b(K(i, j)) = arma::cx_double(1.0, - dt * V(i, j) / 2) - 4. * r;
		}
	}

	// Add the diagonals to the matrices

	A += arma::diagmat(a);
	B += arma::diagmat(b);

// Return the matrices A & B

	return std::make_tuple(A, B);
}

// Define Next_Step function

arma::cx_vec Next_Step(arma::sp_cx_mat A, arma::sp_cx_mat B, arma::cx_vec u)
{
	arma::cx_vec b = B * u;
	arma::cx_vec u_2 = arma::spsolve(A, b);

	return u_2;
}

// Define Initial_State function

arma::cx_vec Initial_State(double h, double x_c, double y_c, double sig_x, double sig_y, double p_x, double p_y)
{
	int M = 1. / h;
	int M_1 = M - 2;
	int M_2 = M_1 * M_1;

	arma::cx_vec u_0(M_2);

	// Create a matrix that translates the indices (i, j) into k

	arma::imat K(M_1, M_1, arma::fill::zeros);

	int k_index = 0;

	for (int i = 0; i < M_1; i++)
	{
		for (int j = 0; j < M_1; j++)
		{
			K(i, j) = k_index;
			k_index += 1;
		}
	}

	// Create the vector u_0

	// Necessary parameters
	std::complex<double> z; // Exponent
	double sig_x2 = sig_x * sig_x;
	double sig_y2 = sig_y * sig_y;

	for (int i = 0; i < M_1; i++)
	{
		for (int j = 0; j < M_1; j++)
		{
			if (i == 0 || j == 0) // Boundrary conditions
			{
				u_0(K(i, j)) = std::complex<double>(0.0, 0.0);
			}
			else
			{
				double x = j * h;
				double y = i * h;
				z = std::complex<double>(-(x - x_c) * (x - x_c) / (2 * sig_x2) - (y - y_c) * (y - y_c) / (2 * sig_y2), p_x * x + p_y * y);
				u_0(K(i, j)) = std::exp(z);
			}
		}
	}

	// Normalise the vector
	double norm = arma::norm(u_0, "fro");
	u_0 = u_0 / norm;

	return u_0;
}

// Define Potential function

// It is important to mention that this function only works for a thickness and aperture of the slits that are multiple of h
// Also the divisions h / thick & h / aperture must result in an even number and the wall that separates the slits must have the same length as the slits

arma::mat Potential(double h, double v_0, int slits, double thick, double aperture)
{
	// Create the Potential matrix without wall

	int M = 1. / h;
	int M_1 = M - 2;
	arma::mat V(M_1, M_1, arma::fill::zeros);
	
	// Create the middle wall

	arma::vec wall(M_1, arma::fill::value(v_0)); // Create a wall with every element v_0 & then we will add the slits as 0
	
	// The idea of the following code is to select the indexes from the middle of the box that will contain a slit (Case odd slits) or wall (Case even slits)
	// This indexes will be displaced to create the slits as a function of the number of slits

	// Select the indexes
	int dy = aperture / h;
	arma::ivec index_y = arma::conv_to<arma::ivec>::from(arma::regspace((M_1 / 2) - (dy / 2), 1 , (M_1 / 2) + (dy / 2 - 1)));
	arma::ivec index;

	// Create the slits
	if (slits % 2 == 0)
	{ 
		for (int i = 0; i < slits / 2; i++)
		{
			index = index_y - dy * (2 * i + 1);
			for (int i = 0; i < index.size(); ++i) 
			{
				wall(index(i)) = 0.;  
			}
		}
		for (int i = 0; i < slits / 2; i++)
		{
			index = index_y + dy * (2 * i + 1);
			for (int i = 0; i < index.size(); ++i)
			{
				wall(index(i)) = 0.;
			}
		}
	}
	if (slits % 2 != 0)
	{
		for (int i = 0; i < (slits + 1) / 2; i++)
		{
			index = index_y - dy * (2 * i);
			for (int i = 0; i < index.size(); ++i)
			{
				wall(index(i)) = 0.;
			}
		}
		for (int i = 0; i < (slits + 1) / 2; i++)
		{
			index = index_y + dy * (2 * i);
			for (int i = 0; i < index.size(); ++i)
			{
				wall(index(i)) = 0.;
			}
		}
	}

	// Add thickness to the wall
	int dx = thick / h;
	arma::ivec index_x = arma::conv_to<arma::ivec>::from(arma::regspace((M_1 / 2) - (dx / 2), 1, (M_1 / 2) + (dx / 2 - 1)));
	
	// Add the wall to the Potential matrix
	for (int i = 0; i < index_x.size(); i++)
	{
		V.col(index_x(i)) = wall;
	}

	return V;
}

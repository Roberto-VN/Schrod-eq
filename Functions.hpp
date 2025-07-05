#include <iostream>
#include <armadillo>
#include <tuple>
#include <complex>

#ifndef __Functions_hpp__   
#define __Functions_hpp__

// Create_Matrix

std::tuple<arma::sp_cx_mat, arma::sp_cx_mat> Create_Matrix(double h, double dt, arma::mat V); // Create the complex matrices A & B that solves the equation

// Next_Step

arma::cx_vec Next_Step(arma::sp_cx_mat A, arma::sp_cx_mat B, arma::cx_vec u); // Update the state of the wave

// Initial_State

arma::cx_vec Initial_State(double h, double x_c, double y_c, double sig_x, double sig_y, double p_x, double p_y); // Create the initial state of the wave

// Potential

arma::mat Potential(double h, double v_0, int slits, double thick, double aperture); // Create the potential matrix

#endif

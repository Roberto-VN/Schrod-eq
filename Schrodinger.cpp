#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <armadillo>
#include <tuple>
#include <complex>

#include "Functions.hpp"

int main(int argc, char* argv[])
{
// Ensure that all parameters are included

    if (argc != 13) {
        std::cerr << "Usage: " << argv[0] << " <dt> <T> <x_c> <sigma_x> <p_x> <y_c> <sigma_y> <p_y> <v_0>\n";
        return 1;
    }

// Parameters set from the terminal

    double h = std::atof(argv[1]);
    double dt = std::atof(argv[2]);  
    double T = std::atof(argv[3]);
    double x_c = std::atof(argv[4]);
    double sig_x = std::atof(argv[5]);
    double p_x = std::atof(argv[6]);
    double y_c = std::atof(argv[7]);
    double sig_y = std::atof(argv[8]);
    double p_y = std::atof(argv[9]);
    double v_0 = std::atof(argv[10]);
    int slits = std::atof(argv[11]);
    int mode = std::atof(argv[12]); // It will define what are we looking for

// Define the dimensions
    int M = 1 / h;
    int M_1 = M - 2;
    int M_2 = M_1 * M_1;

// Create the potential

    // Parameters
    double thick = 0.02;
    double aperture = 0.05;

    arma::mat V = Potential(h, v_0, slits, thick, aperture);

// Set the inital state of the wave

    arma::cx_vec U_0 = Initial_State(h, x_c, y_c, sig_x, sig_y, p_x, p_y);
 
// Create the matrices A & B

    auto[A, B] = Create_Matrix(h, dt, V);

// Run the simulation

    // Set the number of discreticed dt
    int n_dt = T / dt;

    // Temporal vector U
    arma::cx_vec U_n = U_0;

    // Generate the data frame
    std::string filename; 
    int space; // Space for the numbers
    int dec; // Number of decimals

    if (mode == 1) // Meassure the probability deviation along the time
    {
        // Create vector time

        arma::vec time = arma::regspace(0., dt, T);
        
        // Compute the deviation at each time
        
        arma::vec dev(time.size());

        for (int i = 0; i <= n_dt; i++)
        {
            dev(i) = std::abs(1 - arma::norm(U_n, "fro")* arma::norm(U_n, "fro"));
            U_n = Next_Step(A, B, U_n);
        }

        // Set a filename
        filename = "Prob_dev.txt";
            
        // Create the output file
        std::ofstream ofile(filename);
        
        // Necessary parameters for the tables
        space = 27; 
        dec = 16;  

        // Introduce values in the file
        for (int i = 0; i <= n_dt; i++)
        {
            ofile << std::setw(space) << time(i)
                << std::setw(space) << std::setprecision(dec) << dev(i)
                << std::endl;
        }

        // Close the file
        ofile.close();
    }
    if (mode == 2) // Meassure the probabilities of finding the particle at time 0, T/2 & T and also the probability when the particle has been meassured in a certain x position 
    {
        // Create the matrix that will contain the wave values
        arma::cx_mat Wave(3, U_n.size());
        Wave.row(0) = U_n.t();

        // Create the matrices that will contain the wave values & the probabilities
        arma::mat Prob(3, U_n.size());

        for (int i = 0; i < U_n.size(); i++) 
        {
            Prob(0, i) = std::norm(U_n(i));
        }
        for (int k = 0; k < n_dt / 2; k++)
        {
            U_n = Next_Step(A, B, U_n);
        }
        Wave.row(1) = U_n.t();
        for (int i = 0; i < U_n.size(); i++)
        {
            Prob(1, i) = std::norm(U_n(i));
        }
        for (int k = 0; k < n_dt / 2; k++)
        {
            U_n = Next_Step(A, B, U_n);
        }
        Wave.row(2) = U_n.t();
        for (int i = 0; i < U_n.size(); i++)
        {
            Prob(2, i) = std::norm(U_n(i));
        }

    // Probability map

        // Set a filename
        filename = "Prob_map.txt";

        // Create the output file
        std::ofstream ofile(filename);

        // Necessary parameters for the tables
        space = 16;
        dec = 5;

        // Introduce values in the file
        for (int j = 0; j < U_n.size(); j++) {
            for (int i = 0; i < 3; i++) {
                ofile << Prob(i, j) << "\t"; 
            }
            ofile << "\n"; 
        }
            
        // Close the file
        ofile.close();

    // Wave values

        // Set a filename
        filename = "Wave_val.txt";

        // Create the output file
        std::ofstream ofile2(filename);

        // Introduce values in the file
        for (int j = 0; j < U_n.size(); j++) {
            for (int i = 0; i < 3; i++) {
                ofile2 << Wave(i, j) << "\t"; 
            }
            ofile2 << "\n"; 
        }

        // Close the file
        ofile2.close();

    // Probability when the particle has been meassured in a certain x position

        // Create the vector that will contain the probabilities
        arma::vec Prob_y(M_1);

        // Fill the vector with the probabilities 
        double x_meass = 0.8;
        int x_ind = 0.8 / h;

        for (int i = 0; i < M_1; i++)
        {
            Prob_y(i) = Prob(2, x_ind + i * M_1);
        }

        // Normalise the vector
        double norm = arma::norm(Prob_y);
        Prob_y = Prob_y / norm;

    // Create the data frame
        
        // Set a filename
        filename = "Prob_y.txt";

        // Create the output file
        std::ofstream ofile3(filename);

        // Introduce values in the file
        for (int i = 0; i < Prob_y.size(); i++) 
        {
            ofile3 << Prob_y(i) << "\n";
        }

        // Close the file
        ofile3.close();
    }
    if (mode == 3) // Gives back the probabilities of finding the wave in each point of the space at each time
    {
        // Set a filename
        filename = "Wave_sim.txt";

        // Create the output file
        std::ofstream ofile(filename);

        for (int i = 0; i <= n_dt; i++)
        {
            ofile << "Time: " << dt * i << "\n";
            for (int j = 0; j < U_n.size(); j++) 
            {
                ofile << std::norm(U_n(j)) << "\n";
            }
            ofile << "\n\n";

            U_n = Next_Step(A, B, U_n);
        }
        ofile.close();
    }
    if (mode == 4) // Gives back the values of the wave in each point of the space at each time
    {
        // Set a filename
        filename = "Schr_eq.txt";

        // Create the output file
        std::ofstream ofile(filename);

        for (int i = 0; i <= n_dt; i++)
        {
            ofile << "Time: " << dt * i << "\n";
            for (int j = 0; j < U_n.size(); j++)
            {
                ofile << U_n(j) << "\n";
            }
            ofile << "\n\n";

            U_n = Next_Step(A, B, U_n);
        }
        ofile.close();
    }
        
    return 0;
}

/*
 * BremAngleFreqResolved.cpp
 *
 *  Created on: Mar 7, 2020
 *      Author: Evan
 */

/* Numerical relativistic/classical bremsstrahlung yield calculation for either frequency or angle resolved yields.
   Depending on which output is desired, sections of the code must be commented out. These sections are labeled so it is clear which to comment out.*/

// All constants and calculations are performed using atomic units.

#include <iostream>
#include <math.h>
#include <fstream>
#include "eff_z_kr.h"
#include "eff_z_xe.h"
#include "eff_z_U.h"
#include "eff_z_Ar.h"
#include "Vector_Operations.h"
#include <complex>

using namespace std;

// Global Constants
double ke = 1;
double me = 1;
double e = 1;
double	c = 137;
double pi = 3.14159265359;

// Initial Conditions
double X_0 = -10; // Initial position of the electron on the x-axis. The y-axis is the impact parameter, which is a variable in our calculations.
double X_DOT_0 = c * sqrt(1-pow((me * pow(c,2.0)/(me * pow(c,2.0)+1000)), 2.0)); /* Relativistic initial velocity. The "+1000" in the denominator of
																				 the sqrt is the initial kinetic energy of the electron.*/
//double X_DOT_0 = sqrt(2 * 1000/me); // Classical initial velocity; The "1000" term in the sqrt is the initial kinetic energy.
double Y_DOT_0 = 0; // Set to zero for electrons incident from the negative x-axis.
double delta_t =  10 * pow(10, -5); //Initial time step values; This is manipulated in the main() function
double t_start = 0;
double t_final = .45; // Final time of solving differential equation, approximately set to a 1000 Hartree electron traveling a distance of 20 Bohr.


// Acceleration in X direction
double F_x(double x, double y, double vx, double vy, double Z){ /* Relativistic acceleration from the Ion potential on the x-axis calculated by
                                        taking the derivative of the relativistic momentum p=m*v*gamma. T, A, and B are to simplify the EOM. */
	double gamma = 1/sqrt(1 - ((pow(vx,2.0) + pow(vy,2.0))/(pow(c,2.0))));
	   double T = gamma*(1 + (pow(gamma,2.0)/pow(c,2.0))*pow(vy,2.0));
	    double A = (pow(gamma,3.0)/pow(c,2.0)) * vx * vy;
	   double B = gamma*(1 + (pow(gamma,2.0)/pow(c,2.0))*pow(vx,2.0));
	   double Fe_x = (ke*Z*pow(e,2.0)/me)*(x/pow(pow(x,2.0) + pow(y,2.0),1.5));
	   double Fe_y = (ke*Z*pow(e,2.0)/me)*(y/pow(pow(x,2.0) + pow(y,2.0),1.5));
	   double accel_x = ((A/T)*Fe_y - Fe_x)/(B - (pow(A,2.0)/T));
	   //double accel_x = -Fe_x; // Classical acceleration
	   return accel_x;
}

// Acceleration in Y direction
double F_y(double x, double y, double vx, double vy, double Z){ /* Relativistic acceleration from the Ion potential on the y-axis calculated by
                                        taking the derivative of the relativistic momentum p=m*v*gamma. T, A, and B are to simplify the EOM. */
	double gamma = 1/sqrt(1 - ((pow(vx,2.0) + pow(vy,2.0))/(pow(c,2.0))));
	   double T = gamma*(1 + (pow(gamma,2.0)/pow(c,2.0))*pow(vy,2.0));
	    double A = (pow(gamma,3.0)/pow(c,2.0)) * vx * vy;
	   double B = gamma*(1 + (pow(gamma,2.0)/pow(c,2.0))*pow(vx,2.0));
	   double Fe_x = (ke*Z*pow(e,2.0)/me)*(x/pow(pow(x,2.0) + pow(y,2.0),1.5));
	   double Fe_y = (ke*Z*pow(e,2.0)/me)*(y/pow(pow(x,2.0) + pow(y,2.0),1.5));
	   double accel_y = ((A/B)*Fe_x - Fe_y)/(T - (pow(A,2.0)/B));
	   //double accel_y = -Fe_y; // Classical acceleration
	   return accel_y;
}


// trajectory solver
double traj(double b, double freq, double theta, double phi){
	int k;
	double Z;
	int arr_size = (t_final-t_start)/delta_t;
	double x_arr[arr_size];
	double y_arr[arr_size];
	double vx_arr[arr_size];
	double vy_arr[arr_size];
	double ax_arr[arr_size];
	double ay_arr[arr_size];
	int index = 0;
	//ofstream outfile;  // This outputs a trajectory plot if desired.
	//outfile.open("Traj_ar_test.dat");

	// initialize parameters
	double x =  X_0;
	double y = b;
	double vx = X_DOT_0;
	double vy = Y_DOT_0;
	double accelx;
	double accely;
	// While loop solves system
	while(index <= arr_size){ // Numerical integration of differential equation using Riemann sums.
		/* get_z_"species" calls on eff_z_<species>.h to find the effective Z as a function of distance from nucleus using ELSPA electron densities.
		 If calculating a bare nucleus collision Z is set to the atomic number. */
		//Z = get_z_ar(sqrt(pow(x,2.0)+pow(y,2.0)));
		//Z = get_z_kr(sqrt(pow(x,2.0)+pow(y,2.0)));
		//Z = get_z_xe(sqrt(pow(x,2.0)+pow(y,2.0)));
		//Z = get_z_U(sqrt(pow(x,2.0)+pow(y,2.0)));
		Z = 18;
		//Z = 36;
		//Z = 54;
		//Z = 92;
		//outfile << x << " " << y << endl; // This prints to the trajectory plot data file when desired.
		accelx = F_x(x,y,vx,vy,Z);
		accely = F_y(x,y,vx,vy,Z);
		ax_arr[index] = accelx;
		ay_arr[index] = accely;
		vx = vx + (accelx * delta_t);
		vy = vy + (accely * delta_t);
		vx_arr[index] = vx;
		vy_arr[index] = vy;
		x = x + (vx * delta_t);
		y = y + (vy * delta_t);
		x_arr[index] = x;
		y_arr[index] = y;
		index = index + 1;
	}
	//outfile.close(); // Closes the trajectory plot file.

	// Frequency and Angle Resolved Spectra
	/* This part of the function uses equation 14.65 from Jackson's Electrodynamics 2nd edition to calculate the angle and freq resolved spectra. */
	index = 0;
	complex<double> comp_dt = delta_t;
	complex<double> comp_vec[3]{0,0,0};
	while (index <= arr_size){
		double r[3]{x_arr[index], y_arr[index], 0};
		double beta[3]{vx_arr[index]/c, vy_arr[index]/c, 0};
		double beta_dot[3]{ax_arr[index]/c,ay_arr[index]/c,0};
		double n[3]{cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)};
		double *difference = add_vec(n, beta, 1);
		double denominator = pow(1 - dot(beta, n), 2.0);
		double *vector_mess = scale(cross(n, cross(difference, beta_dot)), (1/ denominator));
		complex<double> vector_prod[3]{vector_mess[0], vector_mess[1], vector_mess[2]};
		complex<double> exp_term = exp(1i * freq * ((index * delta_t) - (dot(n,r)/c))) * comp_dt;
		complex<double> *comp_vec_step = scale_comp(vector_prod, exp_term);
		comp_vec[0] = comp_vec[0] + comp_vec_step[0];
		comp_vec[1] = comp_vec[1] + comp_vec_step[1];
		comp_vec[2] = comp_vec[2] + comp_vec_step[2];
		index = index + 1;
	}
	double scalar = pow(comp_vec[0].real(),2.0) + pow(comp_vec[0].imag(),2.0) + pow(comp_vec[1].real(),2.0) + pow(comp_vec[1].imag(),2.0) +
			pow(comp_vec[2].real(),2.0) + pow(comp_vec[2].imag(),2.0);
	double output = (pow(e, 2.0)/(4 * pow(pi,2.0) * c)) * scalar;
	return output;

}


// integration over phi
double step1(double b, double freq, double theta){
	double phi = 0;
	double phi_step = .3142 / 2;
	double energy = 0;
	int index = 0;
	while (index <= 20){
		energy = energy + traj(b, freq, theta, phi)*phi_step*sin(phi);
		phi = phi + phi_step;
		index = index + 1;
	}
	return energy;
}

// integration over theta
double energy_per_freq(double b, double freq){
	double theta = 0;
	double theta_step = .6283 / 2;
	double energy = 0;
	int index = 0;
	while (index <= 20){
		energy = energy + step1(b, freq, theta) * theta_step;
		theta = theta + theta_step;
		index = index + 1;
	}
	return energy;
}


// integration over b (For function of photon energy).
/*
double chi(double bmin, double bmax, double bstep, double initial_chi, double freq){
	int array_size = (bmax-bmin)/bstep;
	double b = bmin;
	int index = 0;
	double chi = initial_chi;
	while (index <= array_size){
		chi = chi + energy_per_freq(b, freq)*2*pi*b*bstep;
		b = b + bstep;
		index = index + 1;
	}
	return chi;
}
*/

//integration over b (For function of phi and theta)
double chi(double bmin, double bmax, double bstep, double initial_chi, double freq, double theta, double phi){
	int array_size = (bmax-bmin)/bstep;
	double b = bmin;
	int index = 0;
	double chi = initial_chi;
	while (index <= array_size){
		chi = chi + traj(b, freq, theta, phi)*2*pi*b*bstep;
		b = b + bstep;
		index = index + 1;
	}
	return chi;
}


// integration over freq (photon energy - for use with angle resolved spectra)
double energyperangle(double gamma_min, double gamma_max, double gamma_step, double theta, double phi){
	int i = 0;
	double freq = gamma_min;
	int array_size = (gamma_max - gamma_min)/gamma_step;
	double freqintegrated = 0;
	while(i <= array_size){
		double ene = 0; // ene is the radiation yield per photon energy
		delta_t = 100 * pow(10, -5.0); // The time step is scaled up to increase computation speed. This is justified given the deflection at the debroglie wavelength is relatively small.
		ene = ene + chi(.138, 10, .4931, ene, freq, theta, phi); // 0.138 is the debroglie wavelength for 1000 Hartree electron, set as the lower limit for the b integration.
		delta_t = 1000 * pow(10, -5.0); // The time step is scaled yet again given the negligible deflection between 10 and 100 bohr.
		ene = ene + chi(10, 100, .9, ene, freq, theta, phi);
		freqintegrated = freqintegrated + ene * gamma_step;
		freq = freq + gamma_step;
		i = i + 1;
	}
	return freqintegrated;
}


// Main Function
int main()
{
	// This is the code to get frequency resolved yields.
	/*
	ofstream outfile;
	outfile.open("freq_resolved.dat");
	int i = 0;
	double freq = 0.1;
	while(i<=100){
		double ene = 0; // ene is the radiation yield per photon energy
		delta_t = 100 * pow(10, -5.0); // The time step is scaled up to increase computation speed. This is justified given the deflection at the debroglie wavelength is relatively small.
		ene = ene + chi(.138, 10, .4931, ene, freq);
		delta_t = 1000 * pow(10, -5.0); // The time step is scaled yet again given the negligible deflection between 10 and 100 bohr.
		ene = ene + chi(10, 100, .9, ene, freq);
		outfile << freq << " " << ene << endl;
		freq = freq + 10;
		i = i + 1;
		if (i == 25){
			cout << "25%" << endl;
		}
		if (i == 50){
			cout << "50%" << endl;
		}
		if (i == 75){
			cout << "75%" << endl;
		}
	}
	*/

	// This is the code to get the angle resolved yield
	ofstream outfile;
	outfile.open("angle_resolved_Z18_Er=1000_repulsion.dat");
	int i = 0;
	double theta_step = (2*pi)/40;
	double phi_step = pi/20;
	double theta = 0;
	while (i <= 40){
		double phi = 0;
		int j = 0;
		while (j <= 20){
			outfile << theta * (180/pi) << " " << phi * (180/pi) << " " << energyperangle(0.1, 1000.1, 10, theta, phi) << endl;
			phi = phi + phi_step;
			j = j + 1;
		}
		theta = theta + theta_step;
		i = i + 1;
		if (i == 10){
			cout << "25%" << endl;
		}
		if (i == 20){
			cout << "50%" << endl;
		}
		if (i == 30){
			cout << "75%" << endl;
		}
	}

	return 0;
}

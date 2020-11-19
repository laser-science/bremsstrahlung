#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

//The code uses atomic units

//Define Constants
double Z  = 10;
double c = 137;
double m = 1;
double e = 1;
double hbar = 1;

double photoIon(double E, double E0, double sig0, double ya, double P, double yw, double y0, double y1){ // This function and its arguments are based off the astrophysics photo-ionization fitting paper
	double x = (E/E0) - y0;
	double y = sqrt(x*x + y1*y1);
	double F = (pow(x-1,2.0)+pow(yw,2.0))*pow(y,(0.5*P)-5.5)*pow(1+sqrt(y/ya),-P);
	double sigma = sig0 * F * (1.0/27.9841); // The last factor is the conversion between Mb and a0^2.
	return sigma;
}

double recomb(double gi, double gf, double Eg, double Er, double E0, double sig0, double ya, double P, double yw, double y0, double y1){
	double sigmaRC = photoIon(Er*27.0 /* Factor of 27 converts Hatree to eV */, E0, sig0, ya, P, yw, y0, y1) * (gi/gf) * (Eg*Eg/(2*m*c*c*Er));
	double Erec = Eg * sigmaRC * 1.0; // the factor of 1 is the fluence
	return Erec;
}

void energy_spectrum_plot_mono(){ // Calculates mono energetic recombination yield
	ofstream ene_spec;
	ene_spec.open("TestCaseRecomb_Ne.dat");
	double ene;
	for(ene=1;ene<=999;ene++){
		ene_spec << ene << " " << recomb(1.0, 4.0, ene, ene-0.78,1.247e1,1.583e3,3.935,7.810,6.558e-2,1.520,1.084e-1) << endl;
	}
}

int main(){
	energy_spectrum_plot_mono();
	return 0;
}

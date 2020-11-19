#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

// This code is designed to calculate the Bethe-Heitler bremsstralung yield using recattering Fluence data.

//The code uses atomic units

//Define Constants
double Z  = 10;
double c = 137;
double m = 1;
double e = 1;
double hbar = 1;
int fluenceSize = 144; // Array size of Fluence data

double brem_per_Egamma(double Er, double Eg){
	double E0 = Er + m*c*c;
	double Ef = E0 - Eg;
	double pf = sqrt(pow(Ef,2.0)-pow(m*pow(c,2.0),2.0));
	double p0 = sqrt(pow(E0,2.0)-pow(m*pow(c,2.0),2.0));
	double mu = m * pow(c,2.0);
	double eps0 = 2 * log((E0 + p0)/mu);
	double eps = 2 * log((Ef + pf)/mu);
	double L = 2 * log((E0*Ef + p0*pf - pow(mu,2.0))/(mu*Eg));
	double alpha = pow(e,2.0)/(hbar*c);
	// Beta-Heitler equation
	double BH = alpha * pow(Z,2.0)*pow(pow(e,2.0)/mu, 2.0) * (pf/p0) * (1/Eg) *(
			(4.0/3.0)-(2*Ef*E0*((pow(pf,2.0)+pow(p0,2.0))/pow(p0*pf,2.0)))+
			pow(mu,2.0)*((eps0*Ef/pow(p0,3.0))+(eps*E0/pow(pf,3.0))-(eps*eps0/(p0*pf)))+
			(((8.0/3.0)*(Ef*E0/(pf*p0)))+(pow(Eg,2.0)/pow(p0*pf,3.0))*(pow(E0*Ef,2.0)+pow(pf*p0,2.0)))*L +
			L*(pow(mu,2.0)*Eg/(2*p0*pf))*(((E0*Ef+pow(p0,2.0))/pow(p0,3.0))*eps0 -((E0*Ef+pow(pf,2.0))/pow(pf,3.0))*eps +
			(2*Eg*E0*Ef/pow(pf*p0,2.0))));
	return BH;
}

// Fluence data is manually added to this array. It can be formated with a coma for easy copy and paste using Origin.
double F_r(int index){
	double list[fluenceSize]{
		3.24192E-218,
		1.65935E-160,
		3.10159E-160,
		4.34803E-160,
		5.4181E-160,
		6.32956E-160,
		7.09857E-160,
		7.73988E-160,
		9.42954E-126,
		1.76253E-125,
		2.47084E-125,
		3.07893E-125,
		3.59688E-125,
		4.03388E-125,
		4.39831E-125,
		1.04972E-102,
		1.9621E-102,
		2.75061E-102,
		3.42755E-102,
		4.00414E-102,
		4.49063E-102,
		2.06102E-86,
		3.85237E-86,
		5.40051E-86,
		6.72961E-86,
		7.86169E-86,
		8.81685E-86,
		2.6148E-74,
		4.88748E-74,
		6.85161E-74,
		8.53783E-74,
		9.97411E-74,
		1.11859E-73,
		5.52225E-65,
		1.0322E-64,
		1.44701E-64,
		1.80312E-64,
		2.10645E-64,
		2.36237E-64,
		1.59178E-57,
		2.9753E-57,
		4.17097E-57,
		5.19748E-57,
		6.07182E-57,
		1.93645E-51,
		3.61952E-51,
		5.07409E-51,
		6.32285E-51,
		1.83724E-46,
		3.43403E-46,
		4.81404E-46,
		5.99879E-46,
		2.07014E-42,
		3.86889E-42,
		5.42343E-42,
		6.75802E-42,
		7.89477E-42,
		8.18406E-39,
		1.52904E-38,
		2.14319E-38,
		2.67044E-38,
		1.10431E-35,
		2.0618E-35,
		2.88928E-35,
		3.59735E-33,
		6.69879E-33,
		9.37904E-33,
		1.16799E-32,
		8.53697E-31,
		1.58549E-30,
		2.21789E-30,
		6.7811E-29,
		1.24812E-28,
		1.74065E-28,
		2.1634E-28,
		4.79703E-27,
		8.77746E-27,
		1.22166E-26,
		2.73726E-25,
		5.00966E-25,
		4.83084E-24,
		8.59204E-24,
		1.18404E-23,
		8.8466E-23,
		1.16539E-6,
		1.08275E-6,
		1.00593E-6,
		9.34528E-7,
		8.68163E-7,
		8.06481E-7,
		7.49155E-7,
		6.95877E-7,
		6.46364E-7,
		6.00351E-7,
		5.57592E-7,
		5.17859E-7,
		4.80938E-7,
		4.46631E-7,
		4.14754E-7,
		3.85137E-7,
		3.5762E-7,
		3.32055E-7,
		3.08304E-7,
		2.8624E-7,
		2.65743E-7,
		2.46703E-7,
		2.29017E-7,
		2.12589E-7,
		1.97331E-7,
		1.83159E-7,
		1.69996E-7,
		1.57772E-7,
		1.46213E-7,
		1.35492E-7,
		1.25549E-7,
		1.16326E-7,
		1.07773E-7,
		9.98445E-8,
		9.24952E-8,
		8.568E-8,
		7.93607E-8,
		7.56212E-8,
		7.20417E-8,
		6.86254E-8,
		6.53747E-8,
		6.22897E-8,
		5.93832E-8,
		5.74238E-8,
		5.55376E-8,
		5.37562E-8,
		5.20919E-8,
		5.4818E-8,
		5.31353E-8,
		5.38676E-8,
		6.13381E-8,
		5.75351E-8,
		7.15895E-8,
		8.51102E-8,
		7.69986E-8,
		1.48822E-7,
		1.5255E-7,
		3.51372E-7,
		1.59838E-6,
		3.39015E-6
	};
	return list[index];
}

// Fluence data is manually added to this array. It can be formated with a coma for easy copy and paste using Origin.
double E_r(int index){
	double list[fluenceSize]{
		0.05095,
		0.05451,
		0.05833,
		0.06241,
		0.06678,
		0.07146,
		0.07646,
		0.08181,
		0.08754,
		0.09366,
		0.10022,
		0.10724,
		0.11474,
		0.12277,
		0.13137,
		0.14056,
		0.1504,
		0.16093,
		0.17219,
		0.18425,
		0.19715,
		0.21095,
		0.22571,
		0.24151,
		0.25842,
		0.27651,
		0.29586,
		0.31657,
		0.33874,
		0.36245,
		0.38782,
		0.41496,
		0.44401,
		0.47509,
		0.50835,
		0.54393,
		0.58201,
		0.62275,
		0.66634,
		0.71299,
		0.76289,
		0.8163,
		0.87344,
		0.93458,
		1,
		1.07,
		1.1449,
		1.22504,
		1.3108,
		1.40255,
		1.50073,
		1.60578,
		1.71819,
		1.83846,
		1.96715,
		2.10485,
		2.25219,
		2.40985,
		2.57853,
		2.75903,
		2.95216,
		3.15882,
		3.37993,
		3.61653,
		3.86968,
		4.14056,
		4.4304,
		4.74053,
		5.07237,
		5.42743,
		5.80735,
		6.21387,
		6.64884,
		7.11426,
		7.61226,
		8.14511,
		8.71527,
		9.32534,
		9.97811,
		10.6766,
		11.4239,
		12.2236,
		13.0793,
		13.9948,
		14.9745,
		16.0227,
		17.1443,
		18.3444,
		19.6285,
		21.0025,
		22.4726,
		24.0457,
		25.7289,
		27.5299,
		29.457,
		31.519,
		33.7253,
		36.0861,
		38.6122,
		41.315,
		44.2071,
		47.3015,
		50.6127,
		54.1555,
		57.9464,
		62.0027,
		66.3429,
		70.9869,
		75.9559,
		81.2729,
		86.962,
		93.0493,
		99.5627,
		106.532,
		113.989,
		121.969,
		130.506,
		139.642,
		149.417,
		159.876,
		171.067,
		183.042,
		195.855,
		209.565,
		224.234,
		239.931,
		256.726,
		274.697,
		293.926,
		314.5,
		336.515,
		360.071,
		385.276,
		412.246,
		441.103,
		471.98,
		505.019,
		540.37,
		578.196,
		618.67,
		661.977,
		708.315,
		757.897,
		810.95
	};
	return list[index];
}


double integrate_over_Er_1(double E_gamma, int index_add){
	int index = 1 + index_add;
	double sum = 0;
	while (index <= fluenceSize-1){
		double delta_Er = E_r(index) - E_r(index - 1);
		double delta_sum = F_r(index) * brem_per_Egamma(E_r(index), E_gamma) * delta_Er;
		sum = sum + delta_sum;
		index = index + 1;
	}
	return sum;
}

double integrate_over_Er_2(double E_gamma, int index_add){
	int index = 1 + index_add;
	double sum = 0;
	while (index <= fluenceSize-1){
		double delta_Er = E_r(index) - E_r(index - 1);
		double delta_sum = F_r(index-1) * brem_per_Egamma(E_r(index-1), E_gamma) * delta_Er;
		sum = sum + delta_sum;
		index = index + 1;
	}
	return sum;
}

void energy_spectrum_plot(){
	ofstream ene_spec;
	ene_spec.open("Energy_Spectrum_U+40.dat");
	int index1 = 1;
	int index2 = 0;
	while (index1 <= fluenceSize-1){
		double Egamma1 = E_r(index1) - .0001; // The reason for the inclusion of the "- 0.0001" is to avoid Egamma being equal to E_r, which would make the bremsstrahlung function undefined.
		double integrate1 = integrate_over_Er_1(Egamma1, index2);
		double Egamma2 = E_r(index1-1) - .0001;
		double integrate2 = integrate_over_Er_2(Egamma2, index2);
		ene_spec << ((Egamma1 + Egamma2)/2) << " " << ((Egamma1 + Egamma2)/2) * ((integrate1 + integrate2)/2) << endl; // Because our integration method necessitates unequal step sizes, either the first or last point of the integration is usually omitted.
		index1 = index1 + 1;                                                                                           // Because of this, we take the average of the two integrations, which is the same method Origin utilizes.
		index2 = index2 + 1;
	}
	ene_spec.close();
}

int main(){
	energy_spectrum_plot();
	return 0;
}

/*
 * Vector_Operations.h
 *
 *  Created on: Mar 20, 2020
 *      Author: Evan
 */

#ifndef VECTOR_OPERATIONS_H_
#define VECTOR_OPERATIONS_H_
#include <complex>

using namespace std;

double dot(double a1[3], double a2[3]){
	double a10 = a1[0];
	double a11 = a1[1];
	double a12 = a1[2];
	double a20 = a2[0];
	double a21 = a2[1];
	double a22 = a2[2];
	double product = a10*a20 + a11*a21 + a12*a22;
	return product;
}


double * cross(double a1[3], double a2[3]){
	static double cross[3];
	double a10 = a1[0];
	double a11 = a1[1];
	double a12 = a1[2];
	double a20 = a2[0];
	double a21 = a2[1];
	double a22 = a2[2];
	cross[0] = a11*a22 - (a21*a12);
	cross[1] = a20*a12 - (a10*a22);
	cross[2] = a10*a21 - (a20*a11);
	return cross;
}


double * scale(double a1[3], double k){
	static double scaled[3];
	double a10 = a1[0];
	double a11 = a1[1];
	double a12 = a1[2];
	scaled[0] = k * a10;
	scaled[1] = k * a11;
	scaled[2] = k * a12;
	return scaled;
}

complex<double> * scale_comp(complex<double> a1[3], complex<double> k){
	static complex<double> scaled[3];
	complex<double> a10 = a1[0];
	complex<double> a11 = a1[1];
	complex<double> a12 = a1[2];
	scaled[0] = k * a10;
	scaled[1] = k * a11;
	scaled[2] = k * a12;
	return scaled;
}

double * add_vec(double a1[3], double a2[3], int k){
	static double added[3];
	double a10 = a1[0];
	double a11 = a1[1];
	double a12 = a1[2];
	double a20 = a2[0];
	double a21 = a2[1];
	double a22 = a2[2];
	double sign;
	if (k == 1){
		sign = -1;
	}
	else{
		sign = 1;
	}
	added[0] = a10 + sign * a20;
	added[1] = a11 + sign * a21;
	added[2] = a12 + sign * a22;
	return added;
}

complex<double> * add_vec_comp(complex<double> a1[3], complex<double> a2[3], int k){
	static complex<double> added[3];
	complex<double> a10 = a1[0];
	complex<double> a11 = a1[1];
	complex<double> a12 = a1[2];
	complex<double> a20 = a2[0];
	complex<double> a21 = a2[1];
	complex<double> a22 = a2[2];
	complex<double> sign;
	if (k == 1){
		sign = -1;
	}
	else{
		sign = 1;
	}
	added[0] = a1[0] + sign * a2[0];
	added[1] = a1[1] + sign * a2[1];
	added[2] = a1[2] + sign * a2[2];
	return added;
}
#endif /* VECTOR_OPERATIONS_H_ */

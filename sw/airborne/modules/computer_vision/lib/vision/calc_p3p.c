/*
 * Copyright (c) 2011, Laurent Kneip, ETH Zurich
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ETH Zurich nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * P3p.cpp
 *
 *  Created on: Nov 2, 2010
 *      Author: Laurent Kneip
 * Description: Compute the absolute pose of a camera using three 3D-to-2D correspondences
 *   Reference: A Novel Parametrization of the P3P-Problem for a Direct Computation of
 *              Absolute Camera Position and Orientation
 *
 *       Input: featureVectors: 3x3 matrix with UNITARY feature vectors (each column is a vector)
 *              worldPoints: 3x3 matrix with corresponding 3D world points (each column is a point)
 *              solutions: 3x16 matrix that will contain the solutions
 *                         form: [ 3x1 position(solution1) 3x3 orientation(solution1) 3x1 position(solution2) 3x3 orientation(solution2) ... ]
 *                         the obtained orientation matrices are defined as transforming points from the cam to the world frame
 *      Output: int: 0 if correct execution
 *                  -1 if world points aligned
 */

#include "P3p.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <complex>
//#include <TooN/SVD.h>
//#include <TooN/determinant.h>

//custom linear algebra functions maybe later to paparazzi math
#define MAT33_TRANS(_mat1,_mat2) {     \
    MAT33_ELMT((_mat1),0,0) = MAT33_ELMT((_mat2),0,0);  \
    MAT33_ELMT((_mat1),0,1) = MAT33_ELMT((_mat2),1,0);  \
    MAT33_ELMT((_mat1),0,2) = MAT33_ELMT((_mat2),2,0);  \
    MAT33_ELMT((_mat1),1,0) = MAT33_ELMT((_mat2),0,1);  \
    MAT33_ELMT((_mat1),1,1) = MAT33_ELMT((_mat2),1,1);  \
    MAT33_ELMT((_mat1),1,2) = MAT33_ELMT((_mat2),2,1);  \
    MAT33_ELMT((_mat1),2,0) = MAT33_ELMT((_mat2),0,2);  \
    MAT33_ELMT((_mat1),2,1) = MAT33_ELMT((_mat2),1,2);  \
    MAT33_ELMT((_mat1),2,2) = MAT33_ELMT((_mat2),2,2);  \
  }



/*P3p::P3p() {
}

P3p::~P3p() {
}
*/

struct FloatVect3 worldPoint_0, worldPoint_1, worldPoint_2, featureVector_0, featureVector_1, featureVector_2,
P1, P2, P3, f1, f2, f3, e1, e2, e3, n1, n2, n3, norm_e3, norm_n1, norm_n3, C, temp1, temp2, temp3, temp4

struct FloatMat33 T, N, R


int P3p_computePoses( TooN::Matrix<3,3> featureVectors, TooN::Matrix<3,3> worldPoints, TooN::Matrix<3,16> & solutions )
{
	// Extraction of world points

	/*TooN::Vector<3> P1 = worldPoints.T()[0];
	TooN::Vector<3> P2 = worldPoints.T()[1];
	TooN::Vector<3> P3 = worldPoints.T()[2];*/
	VECT3_COPY(P1, worldPoint_1);
	VECT3_COPY(P2, worldPoint_2);
	VECT3_COPY(P3, worldPoint_3);
	

	// Verification that world points are not colinear

	/*TooN::Vector<3> temp1 = P2 - P1;
	TooN::Vector<3> temp2 = P3 - P1;*/
	
	VECT3_DIFF(temp1, P2, P1); 
	VECT3_DIFF(temp2, P3, P1); 

	/*if(TooN::norm(temp1 ^ temp2) == 0)
		return -1;*/
	
	VECT3_CROSS_PRODUCT(temp3, temp1, temp2)
	
	if(sqrt(VECT3_NORM2(temp3)) == 0)
	  return -1;

	// Extraction of feature vectors

	/*TooN::Vector<3> f1 = featureVectors.T()[0];
	TooN::Vector<3> f2 = featureVectors.T()[1];
	TooN::Vector<3> f3 = featureVectors.T()[2];*/
	
	VECT3_COPY(f1, featureVector_0);
	VECT3_COPY(f2, featureVector_1);
	VECT3_COPY(f3, featureVector_2);
	
	
	// Creation of intermediate camera frame

	/*TooN::Vector<3> e1 = f1;
	TooN::Vector<3> e3 = f1 ^ f2;
	e3 = e3 / TooN::norm(e3);
	TooN::Vector<3> e2 = e3 ^ e1;*/
	
	VECT3_COPY(e1, f1);
	VECT3_CROSS_PRODUCT(e3, f1, f2);
	norm_e3 = sqrt(VECT3_NORM2(e3));
	VECT3_SDIV(e3, e3, norm_e3);//divide by its norm, to get unit vector. is this possible???

	/*TooN::Matrix<3,3> T;
	T[0] = e1;
	T[1] = e2;
	T[2] = e3;

	f3 = T*f3;*/
	
	MAT33_ROW_VECT3_SMUL(T, 0, e1, 1);
	MAT33_ROW_VECT3_SMUL(T, 1, e2, 1);
	MAT33_ROW_VECT3_SMUL(T, 2, e3, 1);
	
	VECT3_COPY(temp3, f3);//store f3 in temp3 otherwise matrix vector multiplication will get messy?
	
	MAT33_VECT3_MUL(f3, T, temp3); 
	

	// Reinforce that f3[2] > 0 for having theta in [0;pi]

	if( f3.z > 0 )//if( f3[2] > 0 )
	{
		/*f1 = featureVectors.T()[1];
		f2 = featureVectors.T()[0];
		f3 = featureVectors.T()[2];*/
		
		VECT3_COPY(f1, featureVector_1);
		VECT3_COPY(f2, featureVector_0);
		VECT3_COPY(f3, featureVector_2);

		/*e1 = f1;
		e3 = f1 ^ f2;
		e3 = e3 / TooN::norm(e3);
		e2 = e3 ^ e1;*/
		
		VECT3_COPY(e1,f1);
		VECT3_CROSS_PRODUCT(e3, f1, f2);
		norm_e3 = sqrt(VECT3_NORM2(e3));
		VECT3_SDIV(e3, e3, norm_e3);
		VECT3_CROSS_PRODUCT(e2, e3, e1);
		

		/*T[0] = e1;
		T[1] = e2;
		T[2] = e3;*/
		
		
		MAT33_ROW_VECT3_SMUL(T, 0, e1, 1);
		MAT33_ROW_VECT3_SMUL(T, 1, e2, 1);
		MAT33_ROW_VECT3_SMUL(T, 2, e3, 1);

		//f3 = T*f3;
		MAT33_VECT3_MUL(temp3,T,f3);
		VECT3_COPY(f3,temp3);
		

		/*P1 = worldPoints.T()[1];
		P2 = worldPoints.T()[0];
		P3 = worldPoints.T()[2];*/
		
		//TODO make matrix column acessor 
		P1.x = MAT33_ELMT(T,0,1) ;
		P1.y = MAT33_ELMT(T,1,1) ;
		P1.z = MAT33_ELMT(T,2,1) ;
		
		P2.x = MAT33_ELMT(T,0,0) ;
		P2.y = MAT33_ELMT(T,1,0) ;
		P2.z = MAT33_ELMT(T,2,0) ;
		
		P3.x = MAT33_ELMT(T,0,2) ;
		P3.y = MAT33_ELMT(T,1,2) ;
		P3.z = MAT33_ELMT(T,2,2) ;
		
	}

	// Creation of intermediate world frame

	/*TooN::Vector<3> n1 = P2-P1;
	n1 = n1 / TooN::norm(n1);
	TooN::Vector<3> n3 = n1 ^ (P3-P1);
	n3 = n3 / TooN::norm(n3);
	TooN::Vector<3> n2 = n3 ^ n1;*/
	
	VECT3_DIFF(n1, P2, P1);
	norm_n1 = sqrt(VECT3_NORM2(n1));
	VECT3_SDIV(n1, n1, norm_n1);
	//temp before cross product n1 ^ (P3-P1)
	VECT3_DIFF(temp1, P3, P1);
	VECT3_CROSS_PRODUCT(n3, n1, temp1);
	norm_n3 = sqrt(VECT3_NORM2(n3));
	VECT3_SDIV(n3, n3, norm_n3);
	VECT3_CROSS_PRODUCT(n2, n3, n1);

	/*TooN::Matrix<3,3> N;
	N[0] = n1;
	N[1] = n2;
	N[2] = n3;*/
	
	MAT33_ROW_VECT3_SMUL(N, 0, n1, 1);
	MAT33_ROW_VECT3_SMUL(N, 1, n2, 1);
	MAT33_ROW_VECT3_SMUL(N, 2, n3, 1);

	// Extraction of known parameters

	
	
	//P3 = N*(P3-P1);
	VECT3_DIFF(temp1, P3, P1);
	MAT33_VECT3_MUL(P3, N, temp1); 

	/*double d_12 = TooN::norm(P2-P1);
	double f_1 = f3[0]/f3[2];
	double f_2 = f3[1]/f3[2];
	double p_1 = P3[0];
	double p_2 = P3[1];*/
	
	VECT3_DIFF(temp1, P2, P1);
	double f_1 = sqrt(VECT3_NORM2(temp1));
	double f_1 = f3.x/f3.z;
	double f_2 = f3.y/f3.z;
	double p_1 = P3.x;
	double p_2 = P3.y;
	

	/*double cos_beta = f1 * f2;
	double b = 1/(1-pow(cos_beta,2)) - 1;

	if (cos_beta < 0)
		b = -sqrt(b);
	else
		b = sqrt(b);*/
	
	double cos_beta = VECT3_DOT_PRODUCT(f1, f2);
	double b = 1/(1-pow(cos_beta,2)) - 1;

	if (cos_beta < 0)
		b = -sqrt(b);
	else
		b = sqrt(b);
	

	// Definition of temporary variables for avoiding multiple computation

	double f_1_pw2 = pow(f_1,2);
	double f_2_pw2 = pow(f_2,2);
	double p_1_pw2 = pow(p_1,2);
	double p_1_pw3 = p_1_pw2 * p_1;
	double p_1_pw4 = p_1_pw3 * p_1;
	double p_2_pw2 = pow(p_2,2);
	double p_2_pw3 = p_2_pw2 * p_2;
	double p_2_pw4 = p_2_pw3 * p_2;
	double d_12_pw2 = pow(d_12,2);
	double b_pw2 = pow(b,2);

	// Computation of factors of 4th degree polynomial

	//TooN::Vector<5> factors;
	double factors[5];

	factors[0] = -f_2_pw2*p_2_pw4
				 -p_2_pw4*f_1_pw2
				 -p_2_pw4;

	factors[1] = 2*p_2_pw3*d_12*b
				 +2*f_2_pw2*p_2_pw3*d_12*b
				 -2*f_2*p_2_pw3*f_1*d_12;

	factors[2] = -f_2_pw2*p_2_pw2*p_1_pw2
			     -f_2_pw2*p_2_pw2*d_12_pw2*b_pw2
				 -f_2_pw2*p_2_pw2*d_12_pw2
				 +f_2_pw2*p_2_pw4
				 +p_2_pw4*f_1_pw2
				 +2*p_1*p_2_pw2*d_12
				 +2*f_1*f_2*p_1*p_2_pw2*d_12*b
				 -p_2_pw2*p_1_pw2*f_1_pw2
				 +2*p_1*p_2_pw2*f_2_pw2*d_12
				 -p_2_pw2*d_12_pw2*b_pw2
				 -2*p_1_pw2*p_2_pw2;

	factors[3] = 2*p_1_pw2*p_2*d_12*b
				 +2*f_2*p_2_pw3*f_1*d_12
				 -2*f_2_pw2*p_2_pw3*d_12*b
				 -2*p_1*p_2*d_12_pw2*b;

	factors[4] = -2*f_2*p_2_pw2*f_1*p_1*d_12*b
				 +f_2_pw2*p_2_pw2*d_12_pw2
				 +2*p_1_pw3*d_12
				 -p_1_pw2*d_12_pw2
				 +f_2_pw2*p_2_pw2*p_1_pw2
				 -p_1_pw4
				 -2*f_2_pw2*p_2_pw2*p_1*d_12
				 +p_2_pw2*f_1_pw2*p_1_pw2
				 +f_2_pw2*p_2_pw2*d_12_pw2*b_pw2;

	// Computation of roots

	//TooN::Vector<4> realRoots;
	double realRoots[4];

	//this->solveQuartic( factors, realRoots );
	solveQuartic( factors, realRoots );

	// Backsubstitution of each solution

	for(int i=0; i<4; i++)
	{
		double cot_alpha = (-f_1*p_1/f_2-realRoots[i]*p_2+d_12*b)/(-f_1*realRoots[i]*p_2/f_2+p_1-d_12);

		double cos_theta = realRoots[i];
		double sin_theta = sqrt(1-pow(realRoots[i],2));
		double sin_alpha = sqrt(1/(pow(cot_alpha,2)+1));
		double cos_alpha = sqrt(1-pow(sin_alpha,2));

		if (cot_alpha < 0)
			cos_alpha = -cos_alpha;

		/*TooN::Vector<3> C = TooN::makeVector(
				d_12*cos_alpha*(sin_alpha*b+cos_alpha),
				cos_theta*d_12*sin_alpha*(sin_alpha*b+cos_alpha),
				sin_theta*d_12*sin_alpha*(sin_alpha*b+cos_alpha));*/
		
		C.x = d_12*cos_alpha*(sin_alpha*b+cos_alpha);
		C.y = cos_theta*d_12*sin_alpha*(sin_alpha*b+cos_alpha);
		C.z = sin_theta*d_12*sin_alpha*(sin_alpha*b+cos_alpha));

		//C = P1 + N.T()*C;
		MAT33_VECT3_TRANSP_MUL(temp1,N,C);
		VECT3_SUM(C,P1,temp1); 

		/*TooN::Matrix<3,3> R;
		R[0] = TooN::makeVector(	-cos_alpha,		-sin_alpha*cos_theta,	-sin_alpha*sin_theta );
		R[1] = TooN::makeVector(	sin_alpha,		-cos_alpha*cos_theta,	-cos_alpha*sin_theta );
		R[2] = TooN::makeVector(	0,				-sin_theta,				cos_theta );
*/

		//Create matrix TODO if bugs double check this also make matrix assignment macro
		MAT33_ELMT(R, 0, 0) = -cos_alpha;//(row,column)
		MAT33_ELMT(R, 0, 1) = -sin_alpha*cos_theta;
		MAT33_ELMT(R, 0, 2) = -sin_alpha*sin_theta;
		
		MAT33_ELMT(R, 1, 0) = sin_alpha;//(row,column)
		MAT33_ELMT(R, 1, 1) = -cos_alpha*cos_theta;
		MAT33_ELMT(R, 1, 2) = -cos_alpha*sin_theta;
		
		MAT33_ELMT(R, 2, 0) = 0;//(row,column)
		MAT33_ELMT(R, 2, 1) = -sin_theta;
		MAT33_ELMT(R, 2, 2) = cos_theta;
		
		
		//R = N.T()*R.T()*T;
		
		MAT33_TRANS(temp1,R);
		

		/*solutions.T()[i*4] = C;
		solutions.T()[i*4+1] = R.T()[0];
		solutions.T()[i*4+2] = R.T()[1];
		solutions.T()[i*4+3] = R.T()[2];*/
		//copy position result only
		memcpy(result+i*4,C,3*sizeof(float));//will this work???
		
	}

	return 0;
}

int P3p::solveQuartic( TooN::Vector<5> factors, TooN::Vector<4> & realRoots  )
{
	double A = factors[0];
	double B = factors[1];
	double C = factors[2];
	double D = factors[3];
	double E = factors[4];

	double A_pw2 = A*A;
	double B_pw2 = B*B;
	double A_pw3 = A_pw2*A;
	double B_pw3 = B_pw2*B;
	double A_pw4 = A_pw3*A;
	double B_pw4 = B_pw3*B;

	double alpha = -3*B_pw2/(8*A_pw2)+C/A;
	double beta = B_pw3/(8*A_pw3)-B*C/(2*A_pw2)+D/A;
	double gamma = -3*B_pw4/(256*A_pw4)+B_pw2*C/(16*A_pw3)-B*D/(4*A_pw2)+E/A;

	double alpha_pw2 = alpha*alpha;
	double alpha_pw3 = alpha_pw2*alpha;

	/*std::complex<double> P (-alpha_pw2/12-gamma,0);
	std::complex<double> Q (-alpha_pw3/108+alpha*gamma/3-pow(beta,2)/8,0);
	std::complex<double> R = -Q/2.0+sqrt(pow(Q,2.0)/4.0+pow(P,3.0)/27.0);

	std::complex<double> U = pow(R,(1.0/3.0));
	std::complex<double> y;*/
	
	double complex P = -alpha_pw2/12-gamma + 0*I;
	double complex Q = -alpha_pw3/108+alpha*gamma/3-pow(beta,2)/8 + 0*I;
	double complex R = -Q/2.0+csqrt(cpow(Q,2.0)/4.0+cpow(P,3.0)/27.0);

	double complex U = pow(R,(1.0/3.0));
	double complex y;

	if (U.real() == 0)
		y = -5.0*alpha/6.0-pow(Q,(1.0/3.0));
	else
		y = -5.0*alpha/6.0-P/(3.0*U)+U;

	std::complex<double> w = sqrt(alpha+2.0*y);

	std::complex<double> temp;

	temp = -B/(4.0*A) + 0.5*(w+sqrt(-(3.0*alpha+2.0*y+2.0*beta/w)));
	realRoots[0] = temp.real();
	temp = -B/(4.0*A) + 0.5*(w-sqrt(-(3.0*alpha+2.0*y+2.0*beta/w)));
	realRoots[1] = temp.real();
	temp = -B/(4.0*A) + 0.5*(-w+sqrt(-(3.0*alpha+2.0*y-2.0*beta/w)));
	realRoots[2] = temp.real();
	temp = -B/(4.0*A) + 0.5*(-w-sqrt(-(3.0*alpha+2.0*y-2.0*beta/w)));
	realRoots[3] = temp.real();

	return 0;
}

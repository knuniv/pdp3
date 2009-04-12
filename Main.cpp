#include<iostream>
#include "field.h"
#include "E_field.h"
#include "H_field.h"
#include "Time.h"
#include "Particles.h"
#include "Fourier.h"
#include <fstream>
#include <math.h>
using namespace std;
int main() 
{
	PML pml1(0.22,.02,10);
	Geometry geom1(1, 1, 257, 257, &pml1);
	Time time1(0,0, 1E-8,1E-11);
	Particles particle1("Electrons", 1.6e-19, 9.1e-31, 1, &geom1);
	E_field e_field1(&geom1);
	H_field h_field1(&geom1, &particle1);
	Fourier four1(0);
	double** ft;
	ft = new double* [1];
	for (int i=0; i<1;i++)
	{
		ft[i]= new double [geom1.n_grid_2];
	}
	for(int k=0;k<(geom1.n_grid_2);k++)
	{
		ft[0][k]=k;


	}
////   four1.fast_sine_transform(ft, geom1.n_grid_2, 0, true);
// //   four1.fast_fourier_transform(ft, geom1.n_grid_2, 0, true);;
//	e_field1.boundary_conditions();
	e_field1.set_efield_0();
//	e_field1.set_fi_on_z();
//	e_field1.poisson_equation(&geom1);
//	geom1.set_epsilon();
////	e_field1.set_sigma();
//	h_field1.initial_h();
//	particle1.set_j_0();
////	h_field1.magnetostatic_equation(&geom1);
//	ofstream out("test");
//	ofstream out2("test2");



//	for(int i=0;i<=600;i++)
//	{
//		/*	  for (int j=34;j<37;j++)
//		for (int k=62;k<65;k++)*/
//	  //{
//		 //particle1.j2[20][256]=sin(0.2*i);
//		 	 particle1.j2[20][64]=1;
////		 e_field1.e2[j][k]=cos(0.2*i);
//	  //}
//		e_field1.calc_field(&h_field1, &time1, &particle1);
//		h_field1.calc_field(&e_field1, &time1);
//
//
//
//		for(int j=0;(j<geom1.n_grid_1-1);j++)
//		{
//			for(int k=0;k<(geom1.n_grid_2-1);k++)
//			{
//				out<<h_field1.h3[j][k]<<" ";
//////		        out2<<e_field1.fi[j][k]<<" ";
//		};
//			//out<<"\n"; 
////			out2<<"\n";
//		};
//		out<<"\n";
////		out2<<"\n";
//	}
////
//out.close();

//Test relativistic motion
	ofstream out("test");
	particle1.set_v_0();
	particle1.set_x_0();
	for (int i=0; i < 1000000; i++)
	{
		particle1.step_v(&e_field1, &h_field1, &time1);
		particle1.half_step_coord(&time1);
		particle1.half_step_coord(&time1);

		if (i%1000 == 0)
			out<<particle1.v2[0]<<" ";
	}
	out.close();

};
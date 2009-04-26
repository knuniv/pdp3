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
	PML pml1(0.0,0.25, 0.02,10.1);
	Geometry geom1(0.5,0.5, 129, 129, &pml1);
	Time time1(0,0, 1E-8,2E-12);
//	Particles particle1(1000, &geom1);
	E_field e_field1(&geom1);
	H_field h_field1(&geom1);
	Fourier four1(0);
	//double** ft;
	current current1(&geom1);
	charge_density ro1(&geom1);
	//ft = new double* [1];
	//for (int i=0; i<1;i++)
	//{
	//	ft[i]= new double [geom1.n_grid_2];
	//}
	//for(int k=0;k<(geom1.n_grid_2);k++)
	//{
	//	ft[0][k]=k;


	//}
	e_field1.boundary_conditions();
	e_field1.set_efield_0();
	e_field1.set_fi_on_z();
//	e_field1.poisson_equation(&geom1, &ro1);
	geom1.set_epsilon();
//	e_field1.set_sigma();
	h_field1.initial_h();
	ofstream out("test");
	ofstream out2("test2");
	for(int i=0;i<=8000;i++)
	{
			std::cout<<"iteration"<<i;
			std::cout<<endl;
			  for (int j=35;j<37;j++)
				  for (int k=0;k<128;k++)
					 {
//						e_field1.e2[j][k]=sin(0.2*i);
						current1.set_j3(0,k,1);
					 }
		e_field1.calc_field(&h_field1, &time1, &current1, &pml1);
		h_field1.calc_field(&e_field1, &time1);


 if (((i%100)==0)&&(i!=0))
	{
		for(int j=0;(j<geom1.n_grid_1-1);j++)
		{
			for(int k=0;k<(geom1.n_grid_2-1);k++)
			{
				out<<h_field1.h2[j][k]<<" ";
				//out<<geom1.sigma[j][k]<<" ";
////		        out2<<e_field1.fi[j][k]<<" ";
		}
			out<<"\n"; 
//			out2<<"\n";
		}
 
//		out<<"\n";
//		out2<<"\n";
	}
}
//
out.close();


};
#include<iostream>
#include "field.h"
#include "E_field.h"
#include "H_field.h"
#include "Time.h"
#include "Particles.h"
#include "Fourier.h"
#include <fstream>
#include <math.h>
#define  pi = 3.14159265;
using namespace std;
int main() 
{
	PML pml1(0.0,0.25, 0.001,1);
	Geometry geom1(64,512, 65, 513, &pml1);
	Time time1(0,0, 5e-8,1e-10);
	E_field e_field1(&geom1);
	H_field h_field1(&geom1);
	Fourier four1(0);
	bool res;
    ofstream out_coord("coords");
	ofstream out_vel("velocities");
	ofstream out_efield("e_field");
	ofstream out_hfield("h_field");

	current current1(&geom1);
	charge_density rho_new(&geom1);
	charge_density rho_old(&geom1);

	e_field1.boundary_conditions();
	e_field1.set_homogeneous_efield(0.0, 0.0, 0.0);
	h_field1.set_homogeneous_h(0.0, 0.0, 0.0);
	e_field1.set_fi_on_z();
	e_field1.poisson_equation(&geom1, &rho_new);
	geom1.set_epsilon();
//	e_field1.set_sigma();

	Particles new_particles("ions", 1.6e5, 1, 1, &geom1), old_particles("ions", 1, 1, 1, &geom1);
	//////////////////////
	
 	//////////////////////
	double new_x1 =1.4;
	double new_x3 = 1.4;
	double old_x1 = 0.6;
	double old_x3 = 2.1;
    new_particles.x1[0]=old_x1;
	new_particles.x3[0]=old_x3;
	new_particles.charge_weighting(&rho_old);
	new_particles.x1[0]=new_x1;
	new_particles.x3[0]=new_x3;
	new_particles.charge_weighting(&rho_new);
	new_particles.j_weighting(&time1,&current1,new_x1,new_x3,old_x1,old_x3);
    res =  continuity_equation(&time1, &geom1, &current1, &rho_old, &rho_new);
 
	//////////////////////////////////


	new_particles.x1[0] = 52.0001;
    new_particles.x3[0] = 12.0001;
	new_particles.v1[0] = 0.1;
	new_particles.v2[0] = 1e6;
	new_particles.v3[0] = 0.1;
	new_particles.charge_weighting(&rho_new);
	res =  continuity_equation(&time1, &geom1, &current1, &rho_old, &rho_new);

	//0. Half step back
     
    while (time1.current_time < time1.end_time)
	{
        //1. Calculate H field
		h_field1.calc_field(&e_field1, &time1);

		//2. Calculate v
		new_particles.step_v(e_field1.get_field(new_particles.x1[0], new_particles.x3[0]),
			                 h_field1.get_field(new_particles.x1[0], new_particles.x3[0]), &time1);

		//3. Calculate x, calculate J
		old_particles.x1[0] = new_particles.x1[0];
		old_particles.x3[0] = new_particles.x3[0];
		new_particles.charge_weighting(&rho_old);  //continuity equation
        new_particles.step_coord(&time1);
		new_particles.j_weighting(&time1,&current1, new_particles.x1[0], new_particles.x3[0], 
			                                        old_particles.x1[0], old_particles.x3[0]);

        //4. Calculate E
        e_field1.calc_field(&h_field1, &time1, &current1, &pml1);
		
        //continuity equation
		new_particles.charge_weighting(&rho_new);  //continuity equation
		res =  continuity_equation(&time1, &geom1, &current1, &rho_old, &rho_new); 

		out_coord<<new_particles.x1[0]<<" "<<new_particles.x3[0]<<" ";
		out_vel<<new_particles.v1[0]<<" "<<new_particles.v2[0]<<" "<<new_particles.v3[0]<<" ";
		if ((((int)(time1.current_time/time1.delta_t)%100)==0))
		{
			for(int j=0;(j<geom1.n_grid_1-1);j++)
			{
				for(int k=0;k<(geom1.n_grid_2-1);k++)
					{
						out_efield<<e_field1.e3[j][k]<<" ";
						out_hfield<<h_field1.h3[j][k]<<" ";
				    }
	    	}
			out_efield<<"\n"; 
			out_hfield<<"\n"; 
		}
		time1.current_time = time1.current_time + time1.delta_t;
		if (!res)
			break;
 	out_coord<<"\n";
	out_vel<<"\n";

	}

	//particle1.velocity_distribution(1E5);

out_efield.close();
out_hfield.close();
out_vel.close();
out_coord.close();
};
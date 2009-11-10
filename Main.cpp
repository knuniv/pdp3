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

	PML pml1(0.15,0.15, 0.0000001, 0.15);
	Geometry geom1(1.28,1.28, 129, 129, &pml1);
	Time time1(0,0,0,1000e-12,1e-12);
	E_field e_field1(&geom1);
	H_field h_field1(&geom1);
	Fourier four1(0);

	bool res = true;
	int  k, i;
    ofstream out_coord("coords");
	ofstream out_vel("velocities");
	ofstream out_efield("e_field");
	ofstream out_hfield("h_field");
	ofstream curr("curr");

	current current1(&geom1);
	charge_density rho_new(&geom1);
	charge_density rho_old(&geom1);


	//////////////////////////////////////
	double** f_t = new double*[2];
	f_t[0] =new double [9];
	for(k=0;k<9;k++)
		f_t[0][k]=1;
	f_t[0][1]=1;
	four1.dct_2(f_t,9,0,false); 
  // four1.fast_cosine_transform(f_t,9,0,false);
//////////////////////////////////////


	e_field1.boundary_conditions();
	e_field1.set_homogeneous_efield(0.0, 0.0, 0.0);
	h_field1.set_homogeneous_h(0.0, 0.0, 0.0);
	e_field1.set_fi_on_z();
	///////////////////////////////////////////
//	////////////////////////////////////////
//	poisson  equetion testing//
	//for(k=0;k<geom1.n_grid_2;k++)
	//	rho_new.set_ro_weighting(20,k,1e-7);

	//e_field1.poisson_equation2(&geom1,&rho_new);
	////e_field1.poisson_equation(&geom1,&rho_new);
	//bool res2 = e_field1.test_poisson_equation(&rho_new);
//////////////////////////////////////////////////////

////////////////////////////////////////////////
	geom1.set_epsilon() ;
//	e_field1.set_sigma();
    int n_species = 2;
	Particles *new_particles = new Particles[1];
	Particles *old_particles = new Particles[1];

	Particles new_electrons("electrons", -1, 1, 100000, &geom1), old_electrons("electrons", -1, 1, 100000, &geom1);
	Particles new_ions("ions", 1, 1836, 100000, &geom1), old_ions("electrons", 1, 1836, 100000, &geom1);
	new_electrons.load_spatial_distribution(1.6e14, 3.2e14);
	old_electrons.load_spatial_distribution(1.6e14, 3.2e14);
	new_electrons.load_velocity_distribution(0.0);
	old_electrons.load_velocity_distribution(0.0);
	new_ions.load_spatial_distribution(1.6e14, 3.2e14);
	old_ions.load_spatial_distribution(1.6e14, 3.2e14);
	new_ions.load_velocity_distribution(0.0);
	old_ions.load_velocity_distribution(0.0);

	//Particles new_electrons("electrons", -1e5, 1, 1, &geom1), old_electrons("electrons", -1e5, 1, 1, &geom1),
	//	      new_positrons("positrons", 1e5, 1, 1, &geom1), old_positrons("positrons", 1e5, 1, 1, &geom1);

	new_particles[0] = new_electrons;
	old_particles[0] = old_electrons;
    new_particles[1] = new_ions;
    old_particles[1] = old_ions;
	//////////////////////

	
 	//////////////////////


    //two particles test
	//new_particles[0].x1[0] = 0.1500000001;
 //   new_particles[0].x3[0] = .550001;
	//new_particles[0].v1[0] = 0.0;
	//new_particles[0].v2[0] = 0.0;
	//new_particles[0].v3[0] = 0.0e6;

	//new_particles[1].x1[0] = 0.1500000001;
 //   new_particles[1].x3[0] = .650001;
	//new_particles[1].v1[0] = 0.0;
	//new_particles[1].v2[0] = 0.0;
	//new_particles[1].v3[0] = 0.0;
    //////////////////////////


	//0. Half step back

	//weight currents and charges before relaxation period
	for (k=0; k<n_species; k++)
	{
		new_particles[k].azimuthal_j_weighting(&time1, &current1);
		new_particles[k].j_weighting(&time1,&current1,&old_particles[k]);
		new_particles[k].charge_weighting(&rho_new);
	}

	//solve Poisson equation
	e_field1.poisson_equation(&geom1, &rho_new);
	bool res1 = e_field1.test_poisson_equation(&rho_new);
	
	//relaxation period
	while (time1.current_time < time1.relaxation_time)
	{
        //1. Calculate H field
		h_field1.calc_field(&e_field1, &time1);

        //2. Calculate E
        e_field1.calc_field(&h_field1, &time1, &current1, &pml1);
		time1.current_time = time1.current_time + time1.delta_t;
		
       
	}
	time1.current_time = 0.0 ;
     
    while (time1.current_time < time1.end_time)
	{
        //1. Calculate H field
		h_field1.calc_field(&e_field1, &time1);

		//2. Calculate v
		current1.reset_j();
		rho_old.reset_rho();
		for (k=0; k<n_species; k++)
		{
			new_particles[k].step_v(&e_field1, &h_field1, &time1);

		//3. Calculate x, calculate J
			for (i = 0; i< new_particles[k].number; i++)
			{
				old_particles[k].x1[i] = new_particles[k].x1[i];
				old_particles[k].x3[i] = new_particles[k].x3[i];
			}
			new_particles[k].charge_weighting(&rho_old);  //continuity equation
			new_particles[k].half_step_coord(&time1);
			new_particles[k].azimuthal_j_weighting(&time1, &current1);
			new_particles[k].half_step_coord(&time1);
			new_particles[k].j_weighting(&time1,&current1,&old_particles[k]);
		}

        //4. Calculate E
       e_field1.calc_field(&h_field1, &time1, &current1, &pml1);
		
        //continuity equation
		rho_new.reset_rho();
		for (k=0; k<n_species; k++)
			new_particles[k].charge_weighting(&rho_new);  //continuity equation
		res =  continuity_equation(&time1, &geom1, &current1, &rho_old, &rho_new); 

		out_coord<<new_particles[0].x1[0]<<" "<<new_particles[0].x3[0]<<" ";
		out_vel<<new_particles[0].v1[0]<<" "<<new_particles[0].v2[0]<<" "<<new_particles[0].v3[0]<<" ";
		
		if ((((int)(time1.current_time/time1.delta_t))%100==0))
		{
			cout<<time1.current_time<<" ";
			for(int j=0;(j<geom1.n_grid_1-1);j++)
			{
				for(int k=0;k<(geom1.n_grid_2-1);k++)
					{
						out_efield<<e_field1.e3[j][k]<<" ";
						out_hfield<<h_field1.h2[j][k]<<" ";
						curr<<e_field1.e1[j][k]<<" ";
						
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
curr.close();

};
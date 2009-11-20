#include<iostream>
#include "field.h"
#include "E_field.h"
#include "H_field.h"
#include "Time.h"
#include "Particles.h"
#include "Fourier.h"
#include "Poisson.h"
#include "Poisson_neumann.h"
#include "particles_list.h"
#include <fstream>
#include <math.h>
#define  pi = 3.14159265;
using namespace std;
int main() 
{

	PML pml1(0.15,0.15, 0.0000001, 0.15);
	Geometry geom1(1.28,1.28, 129, 129, &pml1);
	Time time1(0,0,0,1e-12,1e-12);
	E_field e_field1(&geom1);
	H_field h_field1(&geom1);
	Fourier four1(0);

	bool res = true;
	int  k, i;
    ofstream out_coord("coords");
	ofstream out_vel("velocities");
	ofstream out_vel1("velocities1");
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
	e_field1.set_homogeneous_efield(0.0, 0.0, 1.0e3);
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
    int n_species = 1;
	Particles *new_particles = new Particles[2];
	Particles *old_particles = new Particles[2];
	particles_list p_list(0);
	Particles electrons("electrons", -1, 1, 1e5, &geom1,&p_list);
	Particles ions("ions", 1, 1836, 1e5, &geom1,&p_list);
	p_list.create_coord_arrays();
	electrons.load_spatial_distribution(1.6e14, 3.2e14);
	
	//electrons.load_velocity_distribution(0.0);

	ions.load_spatial_distribution(1.6e14, 3.2e14);

	//ions.load_velocity_distribution(0.0);

	electrons.velocity_distribution(1000);

	for (i = 0; i< 1e5; i++)
	{
		out_vel1<<electrons.v1[i]<<" "<<electrons.v2[i]<<" "<<electrons.v3[i]<<" ";
	}
	
    //two particles test
	//new_particles[0].x1[0] = 0.0066278076171875002;
 //   new_particles[0].x3[0] = 0.95842972203685850;
	//new_particles[0].v1[0] = 0.0;
	//new_particles[0].v2[0] = 0.0;
	//new_particles[0].v3[0] = 0.0e6;

	//new_particles[1].x1[0] = 0.74797790527343755;
 //   new_particles[1].x3[0] = 0.88477999104412552;
	//new_particles[1].v1[0] = 0.0;
	//new_particles[1].v2[0] = 0.0;
	//new_particles[1].v3[0] = 0.0;
    //////////////////////////
	//0. Half step back

	//p_list.azimuthal_j_weighting(&time1, &current1);
	//p_list.j_weighting(&time1,&current1,);
	//p_list.charge_weighting(&rho_new);

	//weight currents and charges before relaxation period
		//solve Poisson equation
	Poisson_neumann poisson1(&geom1);

	poisson1.poisson_solve(&e_field1, &rho_new);
	bool res1 = poisson1.test_poisson_equation(&e_field1, &rho_new);
	
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
		
			p_list.step_v(&e_field1, &h_field1, &time1);

		//3. Calculate x, calculate J
			p_list.copy_coords();
			p_list.charge_weighting(&rho_old);  //continuity equation
			p_list.half_step_coord(&time1);
			p_list.azimuthal_j_weighting(&time1, &current1);
			p_list.half_step_coord(&time1);
			p_list.j_weighting(&time1,&current1);
		

        //4. Calculate E
       e_field1.calc_field(&h_field1, &time1, &current1, &pml1);
		
        //continuity equation
		rho_new.reset_rho();
		
			p_list.charge_weighting(&rho_new);  //continuity equation
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
out_vel1.close();
out_coord.close();
curr.close();

};
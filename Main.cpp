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
#include "Boundary_Maxwell_conditions.h"
#define  pi = 3.14159265;
using namespace std;
int main() 
{

	PML pml1(0.3,0.0, 0.2, 0.0000001, 0.15);
	Geometry geom1(0.4,3.2, 129, 129, &pml1);
	///////////////////////////////////////////
	ofstream pml("pml");
	for(int i=0;i<geom1.n_grid_1;i++)
		for(int k=0;k<geom1.n_grid_2;k++)
		{
			pml<<geom1.sigma[i][k]<<" ";
		}
		pml.close();
	///////////////////////////////////////////
	Time time1(0,0,0,10000e-12,1e-12);
	E_field e_field1(&geom1);
	H_field h_field1(&geom1);
	Fourier four1(0);
	Boundary_Maxwell_conditions maxwell_rad(&e_field1);
	maxwell_rad.specify_initial_field(&geom1,0,0,0);
	bool res = true;
	int  k, i;
    ofstream out_coord("coords");
	ofstream out_vel("velocities");
	ofstream out_vel1("velocities1");
	ofstream out_e1field("e1_field");
	ofstream out_e3field("e3_field");
	ofstream out_hfield("h_field");
	ofstream curr("curr");

	current current1(&geom1);
	charge_density rho_new(&geom1);
	charge_density rho_old(&geom1);


	e_field1.boundary_conditions();
	e_field1.set_homogeneous_efield(0.0, 0.0, 0);
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
	Particles electrons("electrons", -1, 1, 0*1e6, &geom1,&p_list);
	Particles ions("ions", 1, 1836, 0*1e6, &geom1,&p_list);
	p_list.create_coord_arrays();
	electrons.load_spatial_distribution(2e16, 8e16);


	
	//electrons.load_velocity_distribution(0.0);

	ions.load_spatial_distribution(2e16, 8e16);

		for (i=0; i< 10; i++)
		out_coord<<ions.x1[i]<<" "<<ions.x3[i]<<" ";

	//ions.load_velocity_distribution(0.0);
	
	electrons.velocity_distribution(1e4);
	ions.velocity_distribution(1e3);
	for (i = 0; i< electrons.number; i++)
	{
		out_vel1<<electrons.v1[i]<<" "<<electrons.v2[i]<<" "<<electrons.v3[i]<<" ";
	}
	
    
    //////////////////////////
	//0. Half step back

	//p_list.azimuthal_j_weighting(&time1, &current1);
	//p_list.j_weighting(&time1,&current1,);
	p_list.charge_weighting(&rho_new);

	//weight currents and charges before relaxation period
		//solve Poisson equation

	for(int j=0;(j<geom1.n_grid_1-1);j++)
      for(int k=0;k<(geom1.n_grid_2-1);k++)
		curr<<rho_new.get_ro()[j][k]<<" ";
						
	Poisson_neumann poisson1(&geom1);

	poisson1.poisson_solve(&e_field1, &rho_new);
	bool res1 = poisson1.test_poisson_equation(&e_field1, &rho_new);
	
	//relaxation period
	while (time1.current_time < time1.relaxation_time)
	{
        //1. Calculate H field
		h_field1.calc_field(&e_field1, &time1);

        //2. Calculate E
        e_field1.calc_field(&h_field1, &time1, &current1);
		time1.current_time = time1.current_time + time1.delta_t;
		
       
	}
	time1.current_time = 0.0 ;
     
    while (time1.current_time < time1.end_time)
	{
		//radiation  source
		//maxwell_rad.radiation_source(&geom1,0.4,2e9,0,time1.current_time);
		
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
	   maxwell_rad.probe_mode_exitation(&geom1,&current1, 0.5, 3e8, time1.current_time);
       e_field1.calc_field(&h_field1, &time1, &current1);
	   maxwell_rad.probe_mode_exitation(&geom1,&current1, 0.5, 3e8, time1.current_time);
		
        //continuity equation
		rho_new.reset_rho();
		
			p_list.charge_weighting(&rho_new);  //continuity equation
		res =  continuity_equation(&time1, &geom1, &current1, &rho_old, &rho_new); 

//		out_coord<<new_particles[0].x1[0]<<" "<<new_particles[0].x3[0]<<" ";
//		out_vel<<new_particles[0].v1[0]<<" "<<new_particles[0].v2[0]<<" "<<new_particles[0].v3[0]<<" ";
		
		if  ((((int)(time1.current_time/time1.delta_t))%100==0))
		{
			cout<<time1.current_time<<" ";
			for(int j=0;(j<geom1.n_grid_1-1);j++)
			{
				for(int k=0;k<(geom1.n_grid_2-1);k++)
					{
						out_e1field<<e_field1.e1[j][k]<<" ";
						out_e3field<<e_field1.e3[j][k]<<" ";
						out_hfield<<h_field1.h2[j][k]<<" ";
						//curr<<rho_new.get_ro()[j][k]<<" ";
						
				    }
	    	}
			out_e1field<<"\n"; 
			out_e3field<<"\n";
			out_hfield<<"\n"; 
	        
		}
		time1.current_time = time1.current_time + time1.delta_t;
		if (!res)
			//break;
			cout<<"Error:"<<time1.current_time<<"! ";
 	out_coord<<"\n";
	out_vel<<"\n";
	}

out_e1field.close();
out_e3field.close();
out_hfield.close();
out_vel.close();
out_vel1.close();
out_coord.close();
curr.close();

};
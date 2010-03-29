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
#include "input_output_class.h"
#define  pi = 3.14159265;
using namespace std;
int main() 
{

	PML pml1(0.0,0.0, 0.0, 0.00001, 1);
	Geometry geom1(0.4,6.4, 129, 2049, &pml1);
	double left_plasma_boundary = geom1.second_size*0.2;

	Time time1(0,0,0,100000e-12,1e-12);
	E_field e_field1(&geom1);
	H_field h_field1(&geom1);
	Fourier four1(0);
	input_output_class out_class;
	Boundary_Maxwell_conditions maxwell_rad(&e_field1);
	maxwell_rad.specify_initial_field(&geom1,0,0,0);
	bool res = true;
	int  k, i;
   
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
	particles_list p_list(0);
	Particles electrons("electrons", -1, 1, 1e6, &geom1,&p_list);
	Particles ions("ions", 1, 1836, 1e6, &geom1,&p_list);
	p_list.create_coord_arrays();

	electrons.load_spatial_distribution(0e16, 2.2e15, left_plasma_boundary);
	ions.load_spatial_distribution(0e16, 2.2e15, left_plasma_boundary);

	electrons.velocity_distribution_v2(1e4);
	ions.velocity_distribution_v2(1e3);
	//ofstream out_vel("velocities");
	//ofstream out_coords("coords");
	//for (i = 0; i< electrons.number; i++)
	//{
	//	out_vel<<electrons.v1[i]<<" "<<electrons.v2[i]<<" "<<electrons.v3[i]<<" ";
	//	out_coords<<electrons.x1[i]<<" "<<electrons.x3[i]<<" ";
	//}
	//out_vel.close();
	//out_coords.close();
	   
    //////////////////////////
	//0. Half step back

	//p_list.azimuthal_j_weighting(&time1, &current1);
	//p_list.j_weighting(&time1,&current1,);
	p_list.charge_weighting(&rho_new);
	//electrons.charge_weighting(&rho_new);
	//out_class.out_data("rho",rho_new.get_ro(),1,128,2048);

	//weight currents and charges before relaxation period
		//solve Poisson equation

	/*for(int j=0;(j<geom1.n_grid_1-1);j++)
      for(int k=0;k<(geom1.n_grid_2-1);k++)
		curr<<rho_new.get_ro()[j][k]<<" ";*/
						
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
       e_field1.calc_field(&h_field1, &time1, &current1, &pml1);
		
        //continuity equation
		rho_new.reset_rho();
		
			p_list.charge_weighting(&rho_new);  //continuity equation
		res =  continuity_equation(&time1, &geom1, &current1, &rho_old, &rho_new); 

//		out_coord<<new_particles[0].x1[0]<<" "<<new_particles[0].x3[0]<<" ";
//		out_vel<<new_particles[0].v1[0]<<" "<<new_particles[0].v2[0]<<" "<<new_particles[0].v3[0]<<" ";
		
		if  ((((int)(time1.current_time/time1.delta_t))%100==0))
		{
			cout<<time1.current_time<<" ";
			
			//out_class.out_data("e1",e_field1.e1,100,128,2048);
			out_class.out_data("rho",rho_new.get_ro(),100,geom1.n_grid_1-1,geom1.n_grid_2-1);
			out_class.out_data("e3",e_field1.e3,100,geom1.n_grid_1-1,geom1.n_grid_2-1);
			out_class.out_data("h2",h_field1.h2,100,geom1.n_grid_1-1,geom1.n_grid_2-1);
		}
		time1.current_time = time1.current_time + time1.delta_t;
		if (!res)
			//break;
			cout<<"Error:"<<time1.current_time<<"! ";

	}


};
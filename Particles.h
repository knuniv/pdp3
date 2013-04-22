#pragma once
#include "Geometry.h"
#include "pdp3_time.h"
#include "charge_density.h"
#include"current.h"
#include "triple.h"
#include"particles_list.h"
#include <fstream>
#include<iostream>
using namespace std;
class particles_list;
class E_field;
class H_field;

class Particles
{
public:
	Particles(void);
	Particles(Particles& cp_particles);
	Particles(char* p_name, flcuda p_charge, flcuda p_mass, int p_number,
			  Geometry* geom, particles_list* p_list);
	Particles(char* p_name, double* params,
			  Geometry* geom, particles_list* p_list);
	 ~Particles();
public:
	// The specie name
	char* name;
	// The specie charge
	flcuda charge;
	// The specie mass
	flcuda mass;
	flcuda init_const_mass;
	// Number of particles 
	int number;
	// Particles' coordinates
	flcuda* x1;  //r  
	flcuda* x3;  //z
	// Particles' velocity
	flcuda* v1; //vr
	flcuda* v2; //vphi
	flcuda* v3; //vz

	//indicator if particle is still alive
	int* is_alive;

	flcuda c_light ;
	flcuda c2 ;
    
	//current density
	//temporaty member of Particle class
	//created for the purpose of integration 
	//with exsisting routine of Maxwell equations integration 
	Geometry* geom1;
	particles_list* p_list;
public:
	void charge_weighting(charge_density* ro1);
	void step_v(E_field *e_fld, H_field *h_fld, Time* t);
	virtual void half_step_coord(Time* t);
	void set_j_0();
	void set_v_0();
	void set_x_0();
	flcuda get_gamma(int i);
	flcuda get_gamma_inv(int i);
	void velocity_distribution(flcuda therm_vel);
	void velocity_distribution_v2 (flcuda therm_vel);
	void load_spatial_distribution(double n1, double n2, double left_plasma_boundary,int type);
	void load_velocity_distribution(flcuda v_thermal);
	void simple_j_weighting(Time* time1, current *j1, flcuda x1_new,flcuda x3_new, flcuda x1_old, flcuda x3_old, int i_n, int k_n);
	void simple_constrho_j_weighting(Time* time1, current *j1, flcuda x1_new,flcuda x3_new, flcuda x1_old, flcuda x3_old, int i_n, int k_n);
	void j_weighting(Time* time1, current *j1,flcuda* x1,flcuda* x3);
	void strict_motion_weighting(Time* time1, current *j1, flcuda x1_new,flcuda x3_new, flcuda x1_old, flcuda x3_old);
	void azimuthal_j_weighting(Time* time1, current *j1);
	flcuda  random_reverse(flcuda vel, int power);
	void set_simple_cell(int** cell_arr_jr, int** cell_arr_jz,int start_number, int i_new, int k_new);
	void get_cell_numbers_jr_2(flcuda x1_new,flcuda x3_new, flcuda x1_old, flcuda x3_old,int **cell_arr_jr, int **cell_arr_jz, int* number);
	void get_cell_numbers_jr_1(flcuda x1_new,flcuda x3_new, flcuda x1_old, flcuda x3_old,int *i_return, int* k_return,int* accur);

	};
bool continuity_equation(Time *input_time, Geometry *input_geometry, current *input_J, charge_density *rho_start, charge_density *rho_end);

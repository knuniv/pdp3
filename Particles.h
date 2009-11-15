#pragma once
#include "Geometry.h"
#include "Time.h"
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
	Particles(char* p_name, double p_charge, double p_mass, int p_number,
			  Geometry* geom, particles_list* p_list);
	~Particles();
public:
	// The specie name
	char* name;
	// The specie charge
	double charge;
	// The specie mass
	double mass;
	// Number of particles 
	int number;
	// Particles' coordinates
	double* x1;  //r  
	double* x3;  //z
	// Particles' velocity
	double* v1; //vr
	double* v2; //vphi
	double* v3; //vz

	//indicator if particle is still alive
	bool* is_alive;

	double c_light ;
	double c2 ;
    
	//current density
	//temporaty member of Particle class
	//created for the purpose of integration 
	//with exsisting routine of Maxwell equations integration 
	Geometry* geom1;
	particles_list* p_list;
public:
	void charge_weighting(charge_density* ro1);
	void step_v(E_field *e_fld, H_field *h_fld, Time* t);
	void half_step_coord(Time* t);
	void set_j_0();
	void set_v_0();
	void set_x_0();
	double get_gamma(int i);
	double get_gamma_inv(int i);
	void velocity_distribution(double therm_vel);
	void load_spatial_distribution(double n1, double n2);
	void load_velocity_distribution(double v_thermal);
	void simple_j_weighting(Time* time1, current *j1, double x1_new,double x3_new, double x1_old, double x3_old, int i_n, int k_n);
	void j_weighting(Time* time1, current *j1, Particles *old_part);
	void strict_motion_weighting(Time* time1, current *j1, double x1_new,double x3_new, double x1_old, double x3_old);
	void azimuthal_j_weighting(Time* time1, current *j1);
	//void add_list(particles_list* in_list);
	};
bool continuity_equation(Time *input_time, Geometry *input_geometry, current *input_J, charge_density *rho_start, charge_density *rho_end);

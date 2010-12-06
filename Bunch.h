
#pragma once
#include "Particles.h"
#include "Time.h"

class Bunch :
	public Particles
{
public:
	Bunch(char* p_name, double p_charge, double p_mass, int p_number, Geometry* geom, particles_list* p_list, double b_duration, double b_radius );

	~Bunch(void);
public:
	double duration; //bunch duration
	double radius; //bunch radius
	double n_bunch;// bunch density;
	double vel_bunch; //bunch velocity
public:
	void calc_init_param(double n_b,double b_vel);
	void bunch_inject(Time* time);
	void bunch_inject_calc_E(Geometry* geom, E_field * E_beam,E_field* E, Time* time);
};
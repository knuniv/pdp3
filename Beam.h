#pragma once
#include "Particles.h"
#include "Time.h"

class Beam :
	public Particles
{
public:
	Beam(char* p_name, double p_charge, double p_mass, int p_number, Geometry* geom, particles_list* p_list, double b_duration, double b_radius );

	~Beam(void);
public:
	double duration; //beam duration
	double radius; //beam radius
	double n_beam;// beam density;
	double vel_beam; //beam velocity
public:
	void calc_init_param(double n_b,double b_vel);
	void beam_inject(Time* time);
	void bunch_inject(Time* time, int particles_in_step,double fi,double koef);
	void beam_inject_calc_E(Geometry* geom, E_field * E_beam,E_field* E, Time* time);
};

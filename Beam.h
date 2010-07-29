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
public:
	void calc_init_param(double n_b,double b_vel);
	void beam_inject(double n_b,double b_vel, Time* time);
};

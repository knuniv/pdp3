#pragma once
#include "particles.h"

class Bunch :
	public Particles
{
public:
	Bunch(char* p_name, double p_charge, double p_mass, int p_number, Geometry* geom, particles_list* p_list, double b_radius );

	~Bunch(void);
public:
	double radius; //bunch radius
	double n_bunch;// bunch density;
	double vel_bunch; //bunch velocity
public:
	void calc_init_param(Time* time,int particles_in_step,double n_b,double b_vel);
	void bunch_inject(Time* time, int particles_in_step,double fi,double koef);
	};
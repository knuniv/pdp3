#pragma once
#include"Particles.h"
#include <vector>
#include"charge_density.h"
#include"E_field.h"
#include"H_field.h"
using namespace std;
class Particles;
class particles_list
{
public:
	particles_list(int i);
	~particles_list(void);
public:

public:
	vector<Particles*> part_list;
	void charge_weighting(charge_density* ro1);
	void step_v(E_field *e_fld, H_field *h_fld, Time* t);
	void half_step_coord(Time* t);
	//void set_j_0();
	//void set_v_0();
	//void set_x_0();
	//void velocity_distribution(double therm_vel);
	//void load_spatial_distribution(double n1, double n2);
	//void load_velocity_distribution(double v_thermal);
	void j_weighting(Time* time1, current *j1, Particles *old_part);
	void azimuthal_j_weighting(Time* time1, current *j1);

};

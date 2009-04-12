#pragma once
#include "Geometry.h"
#include "Time.h"
//#include "E_field.h"
//#include "H_field.h"

class E_field;
class H_field;

class Particles
{
public:
	Particles(void);
	Particles(char* p_name, double p_charge, double p_mass, int p_number,
			  Geometry* geom);
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
	double** j1;
	double** j2;
	double** j3;

public:
	double** charge_weighting(Geometry geom);
	void step_v(E_field* e_fld, H_field* H_fld, Time* t);
	void half_step_coord(Time* t);
	void set_j_0();
	void set_v_0();
	void set_x_0();
	double get_gamma(int i);

};

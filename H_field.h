#pragma once
#include "field.h"
#include "Geometry.h"
//#include "E_field.h"
#include "Time.h"
//#include "Particles.h"
#include "Fourier.h"
#include "Triple.h"
class E_field;
class Particles;
class H_field
{
public:
	double** h1;
	double** h2;
	double** h3;
	double** Ar;
	double** Afi;
	double** Az;
	const Particles* particle1;
	const Geometry* geom1;
	H_field(Geometry* geom1, Particles* particle1_t);
	H_field(void);
	~H_field(void);
void calc_field(E_field* e_field1, Time* time1);
void initial_h();
void magnetostatic_equation(Geometry* geom1);
Triple get_field(double x1, double x3);

};

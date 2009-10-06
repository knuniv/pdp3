#pragma once
#include "field.h"
#include "Geometry.h"
//#include "E_field.h"
#include "Time.h"
#include "Particles.h"
#include "Fourier.h"
#include "current.h"
#include "Triple.h"
class E_field;
class H_field
{
public:
	double** h1;
	double** h2;
	double** h3;
	double** h1_half_time;
	double** h2_half_time;
	double** h3_half_time;
	double** Ar;
	double** Afi;
	double** Az;
	const Geometry* geom1;
	H_field(Geometry* geom1);
	H_field(void);
	~H_field(void);
void calc_field(E_field* e_field1, Time* time1);
void set_homogeneous_h(double H1, double H2, double H3);
void magnetostatic_equation(Geometry* geom1);
Triple get_field(double x1, double x3);

};

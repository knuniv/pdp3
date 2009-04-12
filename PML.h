#pragma once
#include"Geometry.h"
class Geometry;
class PML
{
public:
	double comparative_l;
	double sigma1;
	double sigma2;

	void calc_sigma(Geometry* geom1);
	PML(double comp_l_t, double sigma1_t, double sigma2_t);
	PML(void);
	~PML(void);
};

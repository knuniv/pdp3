#pragma once
#include"Geometry.h"
class current
{
public:
	Geometry* geom1;
	current(void);
	current(Geometry* geom1);
	~current(void);
	double** get_j1() const;
	double** get_j2() const;
	double** get_j3() const;
	void set_j1(int i, int k, double value);
	void set_j2(int i, int k, double value);
	void set_j3(int i, int k, double value);
protected:
	double** j1;
	double** j2;
	double** j3;
};

#pragma once
#include"Geometry.h"
class current
{
public:
	Geometry* geom1;
	current(void);
	current(Geometry* geom1);
	~current(void);
	flcuda** get_j1() const;
	flcuda** get_j2() const;
	flcuda** get_j3() const;
	void set_j1(int i, int k, flcuda value);
	void set_j2(int i, int k, flcuda value);
	void set_j3(int i, int k, flcuda value);
	void reset_j();
protected:
	flcuda** j1;
	flcuda** j2;
	flcuda** j3;
};

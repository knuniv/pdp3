#pragma once
#include"E_field.h"
#include"H_field.h"
#include"Geometry.h"

class Boundary_conditions
{
public:
	Boundary_conditions(E_field* ef ,H_field *hf, Geometry *geom);
	~Boundary_conditions(void);
	const Geometry *geom;
	const E_field * ef;
	const H_field *hf;
	void specify_boundary_conditions(double E_fi_upper, double E_fi_left, double E_fi_right, int condition_type);
};

#include "Boundary_conditions.h"
#include"H_field.h"

Boundary_conditions::Boundary_conditions(E_field* ef_t,H_field *hf_t,Geometry *geom_t): ef(ef_t),hf(hf_t),geom(geom_t)
{
}

Boundary_conditions::~Boundary_conditions(void)
{
}

void Boundary_conditions::specify_boundary_conditions(double E_fi_upper, double E_fi_left, double E_fi_right, int condition_type)
{
	int i_max = geom->n_grid_1;
	int k_max = geom->n_grid_2;
	int i,k;
	for (k=0;k<k_max;k++)
	{
		ef->e2[i_max-1][k]=E_fi_upper;
	}
	for (i=0;i<i_max;i++)
	{
		ef->e2[i][0] = E_fi_left;
		ef->e2[i][k_max-1]=E_fi_right;
	}
}
#include "Boundary_Maxwell_conditions.h"

Boundary_Maxwell_conditions::Boundary_Maxwell_conditions(E_field* e_fld_t):e_fld(e_fld_t)
{
}
Boundary_Maxwell_conditions::Boundary_Maxwell_conditions(void)
{
}
Boundary_Maxwell_conditions::~Boundary_Maxwell_conditions(void)
{
}
void Boundary_Maxwell_conditions::specify_initial_field(Geometry* cyl_geom, double E_fi_upper, double E_fi_left, double E_fi_right)
{
int i=0, k=0;
int n_grid1=cyl_geom->n_grid_1;
int n_grid2 = cyl_geom->n_grid_2;
// setazimuthal component electric field initial value
/////////////////////////////////////////////
for (i=0;i<(n_grid1);i++)
{
	e_fld->e2[i][0]=E_fi_left;
	e_fld->e2[i][n_grid2-1]=E_fi_right;
}
	for(k=0;k<n_grid2;k++)
	{
		e_fld->e2[n_grid1-1][k]=E_fi_upper;
	}
/////////////////////////////////////////////

}
void Boundary_Maxwell_conditions::radiation_source(Geometry* cyl_geom,double region_size, double frequency, int wave_type, double time)
{
	int i_max =0, i=0;
	i_max = (int) (region_size/cyl_geom->dr);
	double pi = 3.1415926;

	switch (wave_type)
	{
		// E electromagnetic wave 
		case 0:
    		   {
				for (i=0;i<i_max;i++)
					e_fld->e1[i][0]=e_fld->e1[i][0]*cos(2*pi*frequency*time);
			   }
	   break;
	   // H electromagnetic wave 
		case 1:
    		   {
				for (i=0;i<i_max;i++)
					e_fld->e2[i][0]=e_fld->e2[i][0]*cos(2*pi*frequency*time);
			   }
	   break;
	}
}
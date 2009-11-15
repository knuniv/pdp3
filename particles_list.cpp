#include "particles_list.h"

particles_list::particles_list(int i=0)
{
}

particles_list::~particles_list(void)
{
}
//sicle for all particles kinds in system
void particles_list::charge_weighting(charge_density* rho)
{
	int i=0;
	for(i=0;i<part_list.size();i++)
	{
	 part_list[i]->charge_weighting(rho);
	}
}
void particles_list::step_v(E_field *e_fld, H_field *h_fld, Time* t)
{
	int i=0;
	for(i=0;i<part_list.size();i++)
	{
	 part_list[i]->step_v(e_fld,h_fld,t);
	}
}
void particles_list::half_step_coord(Time* t)
{
	int i=0;
	for(i=0;i<part_list.size();i++)
	{
	 part_list[i]->half_step_coord(t);
	}
}
void particles_list::j_weighting(Time* time1, current *j1, Particles *old_part)
{
	int i = 0;
	for(i=0;i<part_list.size();i++)
	{
		part_list[i]->j_weighting(time1,j1,old_part);
	}
}
void particles_list::azimuthal_j_weighting(Time* time1, current *j1)
{
	int i = 0;
	for(i=0;i<part_list.size();i++)
	{
		part_list[i]->azimuthal_j_weighting(time1,j1);
	}
}
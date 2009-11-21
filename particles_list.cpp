#include "particles_list.h"

particles_list::particles_list(int i=0)
{
	x1_old = 0;
	x3_old = 0;
}

particles_list::~particles_list(void)
{
	for (int i=0; i<(part_list.size());i++)
	{
		delete[] x1_old[i];
		delete[] x3_old[i];
	}
	delete[] x1_old;
	delete[] x3_old;

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
void particles_list::j_weighting(Time* time1, current *j1)
{
	int i = 0;
	for(i=0;i<part_list.size();i++)
	{
		part_list[i]->j_weighting(time1,j1,x1_old[i],x3_old[i]);
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
void particles_list::create_coord_arrays(void)
{
	//creates arrays for storing old particles coordinates
	int kinds_number = part_list.size();
	x1_old = new double* [kinds_number];
	for(int k=0;k<kinds_number;k++)
	{
	x1_old[k]= new double[part_list[k]->number];
	}
	

	x3_old = new double* [kinds_number];
	for(int k=0;k<kinds_number;k++)
	{
	x3_old[k]= new double[part_list[k]->number];
	}
	for(int k=0;k<kinds_number;k++)
		for(int i=0;i<part_list[k]->number;i++)
		{
			x1_old[k][i]=0;
			x3_old[k][i]=0;
		}

}
void particles_list:: copy_coords()
{
	//fuction for copying particles coordinates
	int i=0;
	int k=0;
	int kinds_number = part_list.size();
	for (k=0;k<kinds_number;k++)
		for(i=0;i<part_list[k]->number;i++)
		{
			x1_old[k][i]=part_list[k]->x1[i];
			x3_old[k][i]=part_list[k]->x3[i];
		}
}
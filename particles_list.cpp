#include "particles_list.h"
#include "system_host.cuh"

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
   /* #ifdef BUILD_CUDA
	  TransferEHToCUDA (e_fld->get_1d_e1(), e_fld->get_1d_e2(), e_fld->get_1d_e3(), h_fld->get_1d_h1(), 
		                h_fld->get_1d_h2(), h_fld->get_1d_h3(), e_fld->geom1->n_grid_1, e_fld->geom1->n_grid_2);

    #endif*/
	for(i=0;i<part_list.size();i++)
	{
	 /* #ifdef BUILD_CUDA
		Particles_struct specie = CreateParticles_struct(part_list[i]->charge, part_list[i]->mass, part_list[i]->number, 
			                        e_fld->geom1->n_grid_1, e_fld->geom1->n_grid_2, e_fld->geom1->dr, e_fld->geom1->dz);
		CopySpecie2Cuda (specie);
        TransferXVToCUDA(part_list[i]->x1, part_list[i]->x3, part_list[i]->v1, part_list[i]->v2, part_list[i]->v3, part_list[i]->is_alive, part_list[i]->number);
	    CUDA_StepV(part_list[i]->number, t->delta_t);
	    TransferXVFromCUDA(part_list[i]->x1, part_list[i]->x3, part_list[i]->v1, part_list[i]->v2, part_list[i]->v3, part_list[i]->number);
      #else*/
        part_list[i]->step_v(e_fld,h_fld,t);
     /* #endif*/
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
	x1_old = new flcuda* [kinds_number];
	for(int k=0;k<kinds_number;k++)
	{
	x1_old[k]= new flcuda[part_list[k]->number];
	}
	

	x3_old = new flcuda* [kinds_number];
	for(int k=0;k<kinds_number;k++)
	{
	x3_old[k]= new flcuda[part_list[k]->number];
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
#include "Bunch.h"

Bunch::Bunch(char* p_name, double p_charge, double p_mass, int p_number, Geometry* geom, particles_list* p_list, double b_radius):Particles(p_name,p_charge,p_mass,p_number,geom,p_list)
{
	radius = b_radius;
}
void Bunch::calc_init_param(Time* time,int particles_in_step,double n_b, double b_vel)
{
	n_bunch = n_b;
	vel_bunch = b_vel;
	double dl = vel_bunch*time->delta_t;
	double pi = 3.1415926;
	double n_in_big = 3.1415*radius*radius*dl*n_bunch/particles_in_step;
	charge *= n_in_big;
 	mass *= n_in_big;
	for(int i=0;i<number;i++)
	{
		is_alive[i]=false;
		v1[i]=0;
		v2[i]=0;
		v3[i]=0;
	}
}
void Bunch::bunch_inject(Time* time,int particles_in_step, double fi, double koef)
{
	double dl = vel_bunch*time->delta_t;
	double dr = geom1->dr*1.00000001;
	double dz = geom1->dz*1.00000001;
	double pi = 3.1415926;
	double mod_n = (1 - koef/2.0 + koef/2.0*sin(2*pi*fi*time->current_time));
	int n=0;
	int i=0;
	particles_in_step = particles_in_step*mod_n;
 		while((n<particles_in_step)&&(n<number))
		{
			if(!is_alive[i])
			{
				double	rand_r = random_reverse(i,3);		
				double	rand_z = random_reverse(i,5);
				x1[i] = (radius)*sqrt(rand_r) + dr/2.0;
				x3[i] = dl*(rand_z)-dl;
				v3[i] = vel_bunch;
				is_alive[i] = true;
				n=n+1;
			}
		i=i+1;
		}
		for(int i = 0; i <  number; i++)
		if(x3[i]>(geom1->second_size - dz/2.0))
		{
			is_alive[i]=false;
		}

}	
Bunch::~Bunch(void)
{
}

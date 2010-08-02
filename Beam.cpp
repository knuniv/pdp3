#include "Beam.h"
#include "Time.h"

Beam::Beam(char* p_name, double p_charge, double p_mass, int p_number, Geometry* geom, particles_list* p_list, double b_duration, double b_radius):Particles(p_name,p_charge,p_mass,p_number,geom,p_list)
{
	duration = b_duration;
	radius = b_radius;
};
Beam::~Beam(void)
{
}
void Beam::calc_init_param(double n_b, double b_vel)
{
	double b_lenght = duration*b_vel;
	double n_in_big = 3.1415*radius*radius*b_lenght/number;
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
void Beam::beam_inject(double n_b,double b_vel, Time* time)
{
	double dl = b_vel*time->delta_t;
	int step_num =  duration/time->delta_t;
	int particles_in_step =  number/step_num;
	int start_number = time->current_time/time->delta_t*particles_in_step;
	double dr = geom1->dr*1.00000001;
	double dz = geom1->dz*1.00000001;
    if (time->current_time<duration)
		for(int i = 0; i <  particles_in_step; i++)
		{
		double	rand_r = random_reverse(i,3);		
		double	rand_z = random_reverse(i,5);
			x1[i+start_number] = (radius)*sqrt(rand_r) + dr/2.0;
		
			x3[i+start_number] = dl*(rand_z)+dz/2.0;
			v3[i+start_number] = b_vel;
			is_alive[i+start_number] = true;
		}
}
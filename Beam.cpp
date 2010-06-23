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
double Beam::calc_init_param(double n_b, double b_vel)
{
	double b_lenght = duration*b_vel;
	double n_in_big = 3.1415*radius*radius*b_lenght/number;
	charge *= n_in_big;
	mass *= n_in_big;
}
void Beam::beam_inject(double n_b,double b_vel, Time* time)
{
	double dl = b_vel*time->delta_t;
	int step_num = (int) duration/time->delta_t;
    if (time->current_time<duration)
		for(int i = 0; i < step_num; i++)
		{
		double	rand_r = random_reverse(i,3);		
		double	rand_z = random_reverse(i,5);
			x1[i] = (radius)*sqrt(rand_r) + geom1->dr/2.0;
			x3[i] = dl*sqrt(rand_z)+geom1->dz/2.0;
		}
}
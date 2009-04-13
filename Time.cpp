#include "E_field.h"
#include "H_field.h"
#include "Time.h"

Time::Time(void)
{
}

Time::~Time(void)
{
}

Time::Time(double ct, double st, double et, double dt)
{
	start_time=st;
	current_time=ct;
	end_time=et;
	delta_t=dt;
};
#include "E_field.h"
#include "H_field.h"
#include "pdp3_time.h"

Time::Time(void)
{
}

Time::~Time(void)
{
}

Time::Time(flcuda ct, flcuda st, flcuda rt, flcuda et, flcuda dt)
{
	start_time=st;
	relaxation_time = rt;
	current_time=ct;
	end_time=et;
	delta_t=dt;
};
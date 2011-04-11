#pragma once
#include "particles_struct.h"
class Time
{
public:
	flcuda current_time;
	flcuda start_time;
    flcuda relaxation_time;
	flcuda end_time;
	flcuda delta_t;
	Time(void);
	Time(flcuda ct, flcuda st, flcuda rt, flcuda et, flcuda dt);
	~Time(void);
};

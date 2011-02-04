#pragma once
#include <fstream>
#include<iostream>
#include "particles_struct.h"
using namespace std;
class input_output_class
{
public:
	input_output_class(void);
	~input_output_class(void);
	void out_data(char* comp_name,flcuda** out_value,int step_number,int number, int r_step,int z_step);
	void out_coord(char* comp_name,flcuda* coord_r, flcuda* coord_z,int step_number, int number, int particles_number);
};

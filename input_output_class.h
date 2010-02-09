#pragma once
#include <fstream>
#include<iostream>
using namespace std;
class input_output_class
{
public:
	input_output_class(void);
	~input_output_class(void);
	void out_data(char* comp_name,double** out_value,int number, int r_step,int z_step);

};

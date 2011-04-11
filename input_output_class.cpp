#include "input_output_class.h"
#include <string>
input_output_class::input_output_class(void)
{
}

input_output_class::~input_output_class(void)
{
}
void input_output_class::out_data(char* comp_name, double** out_value,int step_number, int number, int r_step,int z_step)
{
	//static variable for defining number of function calling
    //static int count =0;
	//static variable. increase text file name
	int inc_value = 0;
 // 	if (step_number/number)
	//{
	//	i=i+1;
	////char operations
	//char st_name[50];
	//strcpy(st_name,comp_name);
	//char str_int [50]; 
	//itoa(i,str_int,10);
	//char path[50] = "e:/Plazma/pdp3_result/";
	//
	//strcat(path,st_name);
	//ofstream out_val(path,ios::app);
	//}

   inc_value = step_number/number;
		char st_name[50];
	strcpy(st_name,comp_name);
	char str_int [50]; 
	itoa(inc_value,str_int,10);
	strcat(st_name, str_int);
	char path[50] = "e:/Science[Plasma]/pdp3_result/";
	strcat(path,st_name);
	ofstream out_val(path,ios::app);
	// write  values  into file 
	for (int i=0; i<r_step;i++)
	{
	for(int k=0;k<z_step;k++)
	{
		out_val<<out_value[i][k]<<" ";
	}
	}
	out_val.close();
}

/////
void input_output_class::out_coord(char* comp_name, double* coord_r, double* coord_z,int step_number, int number, int particles_number)
{
	//static variable for defining number of function calling
    //static int count =0;
		//static variable. increase text file name
	int inc_value = 0;
  	/*if (step_number/number)
	{
		i=i+1;
	char st_name[50];
	strcpy(st_name,comp_name);
	char str_int [50]; 
	itoa(i,str_int,10);
	char path[50] = "e:/Plazma/pdp3_result/";
	
	strcat(path,st_name);
	ofstream out_val(path,ios::app);
	}*/
	inc_value = step_number/number;
	char st_name[50];
	strcpy(st_name,comp_name);
	char str_int [50]; 
	itoa(inc_value,str_int,10);
	strcat(st_name, str_int);
	char path[50] = "e:/Science[Plasma]/pdp3_result/";
	strcat(path,st_name);
	ofstream out_val(path,ios::app);
	// write  values  into file 
	for (int i=0; i<particles_number;i++)
	
	{
		out_val<<coord_r[i]<<" ";
		out_val<<coord_z[i]<<" ";
	}

	out_val.close();
}
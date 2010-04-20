#include "input_output_class.h"
#include <string>
input_output_class::input_output_class(void)
{
}

input_output_class::~input_output_class(void)
{
}
void input_output_class::out_data(char* comp_name, double** out_value, int number, int r_step,int z_step)
{
	//static variable for defining number of function calling
    static int count =0;
	//static variable. increase text file name
	static int i=0;
  	if (count/number==i)
	{
		i=i+1;
	//char operations
	char st_name[50];
	strcpy(st_name,comp_name);
	char str_int [50]; 
	itoa(i,str_int,10);
	char path[50] = "e:/Plazma/pdp3_result/";
	
	strcat(path,st_name);
	ofstream out_val(path,ios::app);
	}
	//char add_i = static_cast<char>(i);//
//sprintf(str_int,"%d",i);
	 /*char *add_i= new  char [10];
	 memset(add_i,'\0',10); 
	 itoa(i,add_i,30);
	strcat(comp_name, add_i);*/
	//st_name=st_name+str_int;

		char st_name[50];
	strcpy(st_name,comp_name);
	char str_int [50]; 
	itoa(i,str_int,10);
	strcat(st_name, str_int);
	char path[50] = "e:/Plazma/pdp3_result/";
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
	
	count = count+1;
	out_val.close();
}
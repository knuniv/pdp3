#include "current.h"

current::current(void)
{
}


current::current(Geometry* geom1_t): geom1(geom1_t)
{
	//jr//
	///////////////////////////////////
	//n_grid - кількість ребер
	j1 = new double*[geom1->n_grid_1-1];
	for (int i=0; i<(geom1->n_grid_1-1);i++)
	{
		j1[i]= new double[geom1->n_grid_2];
	}
	///////////////////////////////////
	//jfi//
	j2 = new double*[geom1->n_grid_1];

	for (int i=0; i<(geom1->n_grid_1);i++)
	{
		j2[i]= new double[geom1->n_grid_2];
	}
	//////////////////////////////////////
	//jz//
	//////////////////////////////////////
	j3 = new double*[geom1->n_grid_1];

	for (int i=0; i<(geom1->n_grid_1);i++)
	{
		j3[i]= new double[geom1->n_grid_2-1];
	}
	///////////////////////////////////////


   //////////////////////////////////////
	   //initialization//

	 for (int i=0; i<(geom1->n_grid_1-1);i++)
		for (int k=0; k<(geom1->n_grid_2-1);k++)
		{
			j1[i][k]=0;
			j2[i][k]=0;
			j3[i][k]=0;
		}
	 for (int i=0; i<(geom1->n_grid_1-1);i++)
	  {
		  j1[i][geom1->n_grid_2-1]=0;
	  }

	 for(int k=0;k<(geom1->n_grid_2);k++)
			{
				j2[geom1->n_grid_1-1][k]=0;
			}
	 for(int i=0;i<(geom1->n_grid_1);i++)
		{
		 j2[i][geom1->n_grid_2-1]=0;
		}
	 for(int k=0;k<(geom1->n_grid_2-1);k++)
		{
			j3[geom1->n_grid_1-1][k]=0;
		};
	  

}

current::~current(void)
{
	for (int i=0; i<(geom1->n_grid_1-1);i++)
		delete[]j1[i];
    delete[]j1;

	for (int i=0; i<(geom1->n_grid_1);i++)
		delete[]j2[i];
    delete[]j2;

	for (int i=0; i<(geom1->n_grid_1);i++)
		delete[]j3[i];
    delete[]j3;
}

//////////////////////////////////////////

//////////////////////////////////////////

/////////////////////////////////////////
 //functions for getting access to j arrays//
/////////////////////////////////////////

double** current::get_j1() const
{
	return this->j1;
}

double** current::get_j2() const
{
	return this->j2;
}

double** current::get_j3() const
{
	return this->j3;
}
///////////////////////////////////////////

//////////////////////////////////////////
	// functions for changing values of j//
void current::set_j1(int i, int k, double value)
{
	j1[i][k]+=value;
}

void current::set_j2(int i, int k, double value)
{
	j2[i][k]+=value;
}

void current::set_j3(int i, int k, double value)
{
	j3[i][k]+=value;
}
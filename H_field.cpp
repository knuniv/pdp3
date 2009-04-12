#include "H_field.h"
#include "E_field.h"
#include "Particles.h"
#include <cstdlib>
#include <math.h>
	///Constructor//
////////////////////////////////////////////////

H_field::H_field(Geometry* geom1_l, Particles* particle1_t):geom1(geom1_l), particle1(particle1_t)

{
     //Hr//
/////////////////////////////////////////
	h1 = new double* [geom1->n_grid_1];
	for (int i=0;i<(geom1->n_grid_1);i++)
		{
			h1[i] = new double[geom1->n_grid_2-1];
		}
////////////////////////////////////////

	//Hf//
/////////////////////////////////////////
	h2 = new double* [geom1->n_grid_1-1];
	for (int i=0;i<(geom1->n_grid_1-1);i++)
		{
			h2[i] = new double[geom1->n_grid_2-1];
		}
////////////////////////////////////////

		//Hz//
/////////////////////////////////////////
	h3 = new double* [geom1->n_grid_1-1];
	for (int i=0;i<(geom1->n_grid_1-1);i++)
		{
			h3[i] = new double[geom1->n_grid_2];
		}

////////////////////////////////////////
		//Ar//
/////////////////////////////////////////
	Ar = new double* [geom1->n_grid_1];
	for (int i=0;i<(geom1->n_grid_1);i++)
		{
			Ar[i] = new double[geom1->n_grid_2];
		}

////////////////////////////////////////
		//Afi//
/////////////////////////////////////////
	Afi = new double* [geom1->n_grid_1];
	for (int i=0;i<(geom1->n_grid_1);i++)
		{
			Afi[i] = new double[geom1->n_grid_2];
		}

////////////////////////////////////////
		//Az//
/////////////////////////////////////////
	Az = new double* [geom1->n_grid_1];
	for (int i=0;i<(geom1->n_grid_1);i++)
		{
			Az[i] = new double[geom1->n_grid_2];
		}

////////////////////////////////////
	//initialisation of magnetic potential//
//////////////////////////////////
	for(int i=0;i<geom1->n_grid_1;i++)
		for(int k=0;k<geom1->n_grid_2;k++)
		{
			Ar[i][k]=0;
			Afi[i][k]=0;
			Az[i][k]=0;
		}
//////////////////////////////////
}
////////////////////////////////////////

   //Деструктор///
///////////////////////////////////////
H_field::~H_field(void)
{
	for (int i=0; i<(geom1->n_grid_1);i++)
		delete[]h1[i];
    delete[]h1;

	for (int i=0; i<(geom1->n_grid_1-1);i++)
		delete[]h2[i];
    delete[]h2;

	for (int i=0; i<(geom1->n_grid_1-1);i++)
		delete[]h3[i];
    delete[]h3;

	for (int i=0; i<(geom1->n_grid_1-1);i++)
		delete[]Ar[i];
    delete[]Ar;

	for (int i=0; i<(geom1->n_grid_1-1);i++)
		delete[]Afi[i];
    delete[]Afi;

	for (int i=0; i<(geom1->n_grid_1-1);i++)
		delete[]Az[i];
    delete[]Az;

}

void H_field::initial_h()
{
for (int i=0;i<geom1->n_grid_1;i++)
	for (int k=0;(k<geom1->n_grid_2-1);k++)
	{
			//if (h1[i][k]!= NULL)
				h1[i][k]=0.0;
	}

for (int i=0;i<(geom1->n_grid_1-1);i++)
	for (int k=0;(k<geom1->n_grid_2-1);k++)
	{
			//if (h1[i][k]!= NULL)
				h2[i][k]=0.0;
	}
for (int i=0;i<(geom1->n_grid_1-1);i++)
	for (int k=0;(k<geom1->n_grid_2);k++)
	{
			//if (h1[i][k]!= NULL)
				h3[i][k]=0.0;
	}


		
}
//////////////////////////////////////
	//розрахунок поля//

void H_field::calc_field(E_field* e_field1, Time* time1)
{
	double const_magn0 = 1.26E-6;
	int i=0;
	int k=0;
	//Hr - last i value //
///////////////////////////////////////////////
		for(int k=0;k<(geom1->n_grid_2-1);k++)
		{
			i=geom1->n_grid_1-1;
			this->h1[i][k] = this->h1[i][k]+((e_field1->e2[i][k+1]-e_field1->e2[i][k])/geom1->dz)*time1->delta_t/const_magn0;
		}


	for(int i=0;i<(geom1->n_grid_1-1);i++)
		for(int k=0;k<(geom1->n_grid_2-1);k++)
		{
			this->h1[i][k] = this->h1[i][k]+((e_field1->e2[i][k+1]-e_field1->e2[i][k])/geom1->dz)*time1->delta_t/const_magn0;

			this->h2[i][k] = this->h2[i][k]+((e_field1->e3[i+1][k]-e_field1->e3[i][k])/geom1->dr - (e_field1->e1[i][k+1]-e_field1->e1[i][k])/geom1->dz)*time1->delta_t/const_magn0;

			this->h3[i][k] = this->h3[i][k]-((e_field1->e2[i+1][k]+e_field1->e2[i][k])/(2.0*geom1->dr*(i+0.5))+(e_field1->e2[i+1][k]-e_field1->e2[i][k])/geom1->dr)*time1->delta_t/const_magn0;

		}
///////////////////////////////////////////////
}


//////////////////////////////////////////////
 //magnetostatic equation//
//finding initial disrtibution of magnetic fields//
///////////////////////////////////////////////
void H_field::magnetostatic_equation(Geometry* geom1)
{
	const double pi = 3.141592;
	int i=0;
	int k=0;
	double a=0;
	double c=0;
	double b=0;
	double d=0;
	double* alpha = new double [geom1->n_grid_1];
	double* beta = new double [geom1->n_grid_1];
	Fourier* four1=0;

//////////////////////////////////////////////////
	//


////////////////////////////////////////////////
		for ( k=0; k<(geom1->n_grid_2);k++)
			{
  		
			}
////////////////////////////////////////////////

////////////////////////////////////////////////
	//Ar calculation//
////////////////////////////////////////////////


	//call function for cosine transform//
////////////////////////////////////////////////
////////////////////////////////////////////////
	int temp=geom1->n_grid_2;
	
		for (i=0;i<geom1->n_grid_1-1;i++)
		{
			four1->fast_cosine_transform(particle1->j1, temp, i, false);
		}

		//sweep method//
//////////////////////////////////////////////////////////////////

	for( k=0;k<(geom1->n_grid_2);k++)
	{
		b=2.0/(geom1->dr*geom1->dr)+ pow((pi*k/(geom1->dz*geom1->n_grid_2)),2);
			a = 1.0/(geom1->dr*geom1->dr) + (geom1->dr*(1.0-2.0*(1))-1)/(pow((1.0+2.0*(1)*geom1->dr),2)*geom1->dr);
			c = 1.0/(geom1->dr*geom1->dr) + (geom1->dr*(1.0+2.0*(1))+1)/(pow((1.0+2.0*(1)*geom1->dr),2)*geom1->dr);
			Ar[0][k]=0;//particle1->j1[0][k]*k*geom1->dz*geom1->dr/4;
			d=particle1->j1[1][k]+a*Ar[0][k];
			alpha[2]=c/b;
			beta[2]=d/b;
		for ( i=3;i<(geom1->n_grid_1-2);i++)
			  {	
			//ay-1 - by + cy+1= = -d//
				
				a = 1.0/(geom1->dr*geom1->dr) + (geom1->dr*(1.0-2.0*(i-1))-1)/(pow((1.0+2.0*(i-1)*geom1->dr),2)*geom1->dr);
				c = 1.0/(geom1->dr*geom1->dr) + (geom1->dr*(1.0+2.0*(i-1))+1)/(pow((1.0+2.0*(i-1)*geom1->dr),2)*geom1->dr);
				d = particle1->j1[i-1][k];

				alpha[i] = c/(b-alpha[i-1]*a);
				beta[i]  = (d+beta[i-1]*a)/(b-alpha[i-1]*a);

			  }
				i = geom1->n_grid_1-2;
				a = 1.0/(geom1->dr*geom1->dr) + (geom1->dr*(1.0-2.0*(i-1))-1)/(pow((1.0+2.0*(i-1)*geom1->dr),2)*geom1->dr);
				c = 1.0/(geom1->dr*geom1->dr) + (geom1->dr*(1.0+2.0*(i-1))+1)/(pow((1.0+2.0*(i-1)*geom1->dr),2)*geom1->dr);
				d = particle1->j1[i-1][k];

				Ar[geom1->n_grid_1-3][k]=-d-c*Ar[geom1->n_grid_1-2][k]-a*beta[geom1->n_grid_1-3]/(a*alpha[geom1->n_grid_1-3]-b);

				for(i=(geom1->n_grid_1-4);i>=1;i--)
				 {
					Ar[i][k]=beta[i+1]+alpha[i+1]*Ar[i+1][k];
				 }
      }


	 temp=geom1->n_grid_2;
			for (i=0;i<geom1->n_grid_1-1;i++)
		{
			four1->fast_cosine_transform(Ar, temp, i, true);
		}

}

Triple H_field::get_field(double x1, double x3)
{
	Triple components(0.0, 0.0, 0.0);

	return components;
}
#include "H_field.h"
#include "E_field.h"
#include <cstdlib>
#include <math.h>

	///Constructor//
////////////////////////////////////////////////

H_field::H_field(Geometry* geom1_l):geom1(geom1_l)

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

	     //Hr_half_time//
/////////////////////////////////////////
	h1_half_time = new double* [geom1->n_grid_1];
	for (int i=0;i<(geom1->n_grid_1);i++)
		{
			h1_half_time[i] = new double[geom1->n_grid_2-1];
		}
////////////////////////////////////////

	//Hf half time//
/////////////////////////////////////////
	h2_half_time = new double* [geom1->n_grid_1-1];
	for (int i=0;i<(geom1->n_grid_1-1);i++)
		{
			h2_half_time[i] = new double[geom1->n_grid_2-1];
		}
////////////////////////////////////////

		//Hz half time//
/////////////////////////////////////////
	h3_half_time = new double* [geom1->n_grid_1-1];
	for (int i=0;i<(geom1->n_grid_1-1);i++)
		{
			h3_half_time[i] = new double[geom1->n_grid_2];
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
	for(int i=0;i<geom1->n_grid_1-1;i++)
		for(int k=0;k<geom1->n_grid_2-1;k++)
		{
			h1[i][k]=0;
			h2[i][k]=0;
			h3[i][k]=0;
			h1_half_time[i][k]=0;
			h2_half_time[i][k]=0;
			h3_half_time[i][k]=0;
		}

		for(int k=0;k<geom1->n_grid_2-1;k++)
		{
			h1[geom1->n_grid_1-1][k]=0;
			h1_half_time[geom1->n_grid_1-1][k]=0;
		}


//////////////////////////////////
}
////////////////////////////////////////

   //����������///
///////////////////////////////////////
H_field::~H_field(void)
{
	for (int i=0; i<(geom1->n_grid_1);i++)
	{
		delete[]h1[i];
		delete[]h1_half_time[i];
	}
    delete[]h1;
	delete[]h1_half_time;

	for (int i=0; i<(geom1->n_grid_1-1);i++)
	{
		delete[]h2[i];
		delete[]h2_half_time[i];
	}
    delete[]h2;
	delete[]h2_half_time;

	for (int i=0; i<(geom1->n_grid_1-1);i++)
	{
		delete[]h3[i];
		delete[]h3_half_time[i];
	}
    delete[]h3;
	delete[]h3_half_time;

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
	//���������� ����//

void H_field::calc_field(E_field* e_field1, Time* time1)
{
	double const_magn0 = 1.26E-6;
	double alpha;
	int i=0;
	int k=0;
	//Hr - last i value //
///////////////////////////////////////////////
		for(int k=0;k<(geom1->n_grid_2-1);k++)
		{
			i=geom1->n_grid_1-1;
			alpha=((e_field1->e2[i][k+1]-e_field1->e2[i][k])/geom1->dz)/const_magn0;
	
			this->h1_half_time[i][k]=this->h1[i][k]+alpha*time1->delta_t/2;
			this->h1[i][k] = this->h1[i][k]+alpha*time1->delta_t;
		}


	for(int i=0;i<(geom1->n_grid_1-1);i++)
		for(int k=0;k<(geom1->n_grid_2-1);k++)
		{
			alpha=((e_field1->e2[i][k+1]-e_field1->e2[i][k])/geom1->dz)/const_magn0;
	
			this->h1_half_time[i][k]=this->h1[i][k]+alpha*time1->delta_t/2;
			this->h1[i][k] = this->h1[i][k]+alpha*time1->delta_t;

			alpha=((e_field1->e3[i+1][k]-e_field1->e3[i][k])/geom1->dr - (e_field1->e1[i][k+1]-e_field1->e1[i][k])/geom1->dz)/const_magn0;
			this->h2_half_time[i][k] = this->h2[i][k]+alpha*time1->delta_t/2;
			this->h2[i][k] = this->h2[i][k]+alpha*time1->delta_t;

			alpha= ((e_field1->e2[i+1][k]+e_field1->e2[i][k])/(2.0*geom1->dr*(i+0.5))+(e_field1->e2[i+1][k]-e_field1->e2[i][k])/geom1->dr)/const_magn0;
			this->h3_half_time[i][k] = this->h3[i][k]-alpha*time1->delta_t/2;
			this->h3[i][k] = this->h3[i][k]-alpha*time1->delta_t;

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
//	int temp=geom1->n_grid_2;
//	
//		for (i=0;i<geom1->n_grid_1-1;i++)
//		{
//			four1->fast_cosine_transform(particle1->j1, temp, i, false);
//		}
//
//		//sweep method//
////////////////////////////////////////////////////////////////////
//
//	for( k=0;k<(geom1->n_grid_2);k++)
//	{
//		b=2.0/(geom1->dr*geom1->dr)+ pow((pi*k/(geom1->dz*geom1->n_grid_2)),2);
//			a = 1.0/(geom1->dr*geom1->dr) + (geom1->dr*(1.0-2.0*(1))-1)/(pow((1.0+2.0*(1)*geom1->dr),2)*geom1->dr);
//			c = 1.0/(geom1->dr*geom1->dr) + (geom1->dr*(1.0+2.0*(1))+1)/(pow((1.0+2.0*(1)*geom1->dr),2)*geom1->dr);
//			Ar[0][k]=0;//particle1->j1[0][k]*k*geom1->dz*geom1->dr/4;
//			d=particle1->j1[1][k]+a*Ar[0][k];
//			alpha[2]=c/b;
//			beta[2]=d/b;
//		for ( i=3;i<(geom1->n_grid_1-2);i++)
//			  {	
//			//ay-1 - by + cy+1= = -d//
//				
//				a = 1.0/(geom1->dr*geom1->dr) + (geom1->dr*(1.0-2.0*(i-1))-1)/(pow((1.0+2.0*(i-1)*geom1->dr),2)*geom1->dr);
//				c = 1.0/(geom1->dr*geom1->dr) + (geom1->dr*(1.0+2.0*(i-1))+1)/(pow((1.0+2.0*(i-1)*geom1->dr),2)*geom1->dr);
//				d = particle1->j1[i-1][k];
//
//				alpha[i] = c/(b-alpha[i-1]*a);
//				beta[i]  = (d+beta[i-1]*a)/(b-alpha[i-1]*a);
//
//			  }
//				i = geom1->n_grid_1-2;
//				a = 1.0/(geom1->dr*geom1->dr) + (geom1->dr*(1.0-2.0*(i-1))-1)/(pow((1.0+2.0*(i-1)*geom1->dr),2)*geom1->dr);
//				c = 1.0/(geom1->dr*geom1->dr) + (geom1->dr*(1.0+2.0*(i-1))+1)/(pow((1.0+2.0*(i-1)*geom1->dr),2)*geom1->dr);
//				d = particle1->j1[i-1][k];
//
//				Ar[geom1->n_grid_1-3][k]=-d-c*Ar[geom1->n_grid_1-2][k]-a*beta[geom1->n_grid_1-3]/(a*alpha[geom1->n_grid_1-3]-b);
//
//				for(i=(geom1->n_grid_1-4);i>=1;i--)
//				 {
//					Ar[i][k]=beta[i+1]+alpha[i+1]*Ar[i+1][k];
//				 }
//      }
//
//
//	 temp=geom1->n_grid_2;
//			for (i=0;i<geom1->n_grid_1-1;i++)
//		{
//			four1->fast_cosine_transform(Ar, temp, i, true);
//		}

}

////////////////////////////////////////////////////////////
	//function for magnetic field weighting//

Triple H_field::get_field(double x1, double x3, int half)
{

	{	
	int i_r=0;  // number of particle i cell 
	int k_z=0;  // number of particle k cell
	
	double pi = 3.14159;
	double dr = geom1->dr;
	double dz = geom1->dz;
	double r1, r2, r3; // temp variables for calculation
	double dz1, dz2;   // temp var.: width of k and k+1 cell 
	double hr =0;
	double hfi =0;
	double hz =0;
	double vol_1 =0; //  volume of i cell; Q/V, V - volume of elementary cell 
	double vol_2 =0; //  volume of i+1 cell;
	
	double value =0;
////////////////////////
	r1 = x1-0.5*dr;
	r3 = x1+0.5*dr;
///////////////////////

    // weighting of H_z//
///////////////////////////////////////////////////
	//finding number of cell. example dr=0.5,  x1 = 0.7, i_r =0;!!
	 i_r = (int)ceil((x1-0.5*dr)/geom1->dr)-1;
	 k_z = (int)ceil((x3)/geom1->dz)-1;

	 vol_1 = pi*dz*dr*dr*(2*i_r+1);
	 vol_2 = pi*dz*dr*dr*(2*i_r+3);
	 dz1 = (k_z+1)*dz-x3;
	 dz2 = x3 - k_z*dz;
	 r2 = (i_r+1)*dr;
    ///////////////////////////////////////
  
    //weighting Hz[i][k]//
   hz = hz + h3[i_r][k_z]*(pi*dz1*(r2*r2-r1*r1))/vol_1;

	//weighting Hz[i+1][k]//
   hz = hz + h3[i_r+1][k_z]*(pi*dz1*(r3*r3-r2*r2))/vol_2;

   //weighting Hz[i][k+1]//
   hz= hz + h3[i_r][k_z+1]*(pi*dz2*(r3*r3-r2*r2))/vol_1;

   //weighting Hz[i+1][k+1]//
  hz = hz + h3[i_r+1][k_z+1]*(pi*dz2*(r3*r3-r2*r2))/vol_2;
   
///////////////////////////////////////////////////////



	     // weighting of Hr//
///////////////////////////////////////////////////////

  //finding number of cell. example dz=0.5,  x3 = 0.7, z_k =0;!!
	i_r = (int)ceil((x1)/geom1->dr)-1;
	k_z = (int)ceil((x3-0.5*dz)/geom1->dz)-1;
  ////////////////////////////

   if(x1>dr/2)
	{
		vol_1 = pi*dz*dr*dr*2*i_r;
    }
   else 
   {
        vol_1 = pi*dz*dr*dr/4.0; //volume of first cell
   }

		  r2 = (i_r+0.5)*dr;

		  vol_2 = pi*dz*dr*dr*(2*i_r+2);
		  dz1 = (k_z+1.5)*dz - x3;
		  dz2 = x3 - (k_z+0.5)*dz;
		  //////////////////////////////////////

		   //weighting Hr[i][k]//
		   hr = hr + h1[i_r][k_z]*(pi*dz1*(r2*r2-r1*r1))/vol_1;

		  //weighting Hr[i+1][k]//
		   hr = hr + h1[i_r+1][k_z]*pi*dz1*(r3*r3-r2*r2)/vol_2;   

          //weighting Hr[i][k+1]//
		   hr = hr + h1[i_r][k_z+1]*pi*dz2*(r2*r2-r1*r1)/vol_1;
   
         //weighting Hr[i+1][k+1]//
		   hr = hr + h1[i_r+1][k_z+1]*pi*dz2*(r3*r3-r2*r2)/vol_2;    
    
  
///////////////////////////////////////////////////////

	 // weighting of H_fi//
///////////////////////////////////////////////////////
 
	
	     //finding number of cell. example dz=0.5,  x3 = 0.7, z_k =0;
		  i_r = (int)ceil((x1-0.5*dr)/geom1->dr)-1;
		  k_z = (int)ceil((x3-0.5*dz)/geom1->dz)-1;

		  r2 = (i_r+1)*dr;
		  vol_1 = pi*dz*dr*dr*(2*i_r+1);
		  vol_2 = pi*dz*dr*dr*(2*i_r+3);
		  dz1 = (k_z+1.5)*dz-x3;
		  dz2 = x3-(k_z+0.5)*dz;
		  //////////////////////////////////////
		  //weighting Hfi[i][k]//
		  hfi = hfi + h2[i_r][k_z]*pi*dz1*(r2*r2-r1*r1)/vol_1;

		  //weighting Hfi[i+1][k]//
		  hfi = hfi + h2[i_r+1][k_z]*pi*dz1*(r3*r3-r2*r2)/vol_2;
		   
          //weighting Hfi[i][k+1]//
		   hfi = hfi + h2[i_r][k_z+1]*dz2*pi*(r2*r2-r1*r1)/vol_1;
   
         //weighting Hfi[i+1][k+1]//
		   hfi = hfi + h2[i_r+1][k_z+1]*pi*dz2*(r3*r3-r2*r2)/vol_2;
  
		   
	Triple components(hr, hfi, hz);

	return components;
}

}
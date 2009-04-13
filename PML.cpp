#include "PML.h"
#include "math.h"

PML::PML(void)
{
}

PML::PML(double comp_l_t, double sigma1_t, double sigma2_t)
{
	comparative_l = comp_l_t;
	sigma1 = sigma1_t;
	sigma2 = sigma2_t;
		
};

PML::~PML(void)
{
}
		//sigma assignment//
//////////////////////////////////////////////////////////////
void PML::calc_sigma(Geometry* geom1)
{	
	int i=0;
	int k=0;

//////////////////////////////////////////

  // defining lenght of sigma calculation//
double lenght_sigma_z =geom1->dz * (floor(geom1->n_grid_2*comparative_l));
	//sigma on r=0 and r=r axis//
/////////////////////////////////////////////
 for(k=0;k<geom1->n_grid_2;k++)
	 for(i=0;i<geom1->n_grid_1;i++)
	 {
			if (geom1->dz*k<=lenght_sigma_z)
				{
					geom1->sigma[i][k] = sigma1 + (sigma2-sigma1)/(pow(lenght_sigma_z,2))*pow((lenght_sigma_z-geom1->dz*k),2);
				}

			if ((geom1->second_size - geom1->dz*(k))<=lenght_sigma_z)
				{
					 geom1->sigma[i][k] = sigma1 + (sigma2-sigma1)/(pow(lenght_sigma_z,2))*pow((geom1->dz*(k+1)-geom1->second_size+lenght_sigma_z),2);
				}
	 }
//////////////////////////////////////////////	
// defining lenght of sigma calculation//
double lenght_sigma_r =geom1->dr * (floor(geom1->n_grid_1*comparative_l));
    // sigma assigning//	 
	//sigma on z=z axis//
//////////////////////////////////////////////
 for(i=0;i<geom1->n_grid_1;i++)
	 for(k=0;k<geom1->n_grid_2;k++)
  if ((geom1->first_size - geom1->dr*(i)) <= lenght_sigma_r)
//	  if ((geom1->dz*(k+1) > lenght_sigma_r) && (geom1->second_size - geom1->dz*(k+1) > lenght_sigma_z))
//		{
//			 geom1->sigma[i][k] = sigma1 + (sigma2-sigma1)/(pow(lenght_sigma_r,2))*pow((geom1->dr*(i+1)-geom1->first_size+lenght_sigma_r),2);
//		}
	   if (((geom1->first_size - geom1->dr*(i+1)) < lenght_sigma_r * geom1->dz * (k)/lenght_sigma_z)  &&  ( (geom1->first_size - geom1->dr*(i)) <= (lenght_sigma_r/lenght_sigma_z*(geom1->second_size - geom1->dz*(k)))))
				geom1->sigma[i][k] = sigma1 + (sigma2-sigma1)/(pow(lenght_sigma_r,2))*pow((geom1->dr*(i+1)-geom1->first_size+lenght_sigma_r),2);
/////////////////////////////////////////////
return ;
}
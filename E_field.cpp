#include "E_field.h"
#include "H_field.h"
#include "Math.h"
#include "Fourier.h"

E_field::E_field(): epsilon0(8.85E-12)
{
}
////конструктор. ініціалізація динамічних массивів
E_field::E_field(Geometry* geom1_t): geom1(geom1_t), epsilon0(8.85E-12)
{
	//Er//
	///////////////////////////////////
	//n_grid - кількість ребер
	e1 = new double*[geom1->n_grid_1-1];
	for (int i=0; i<(geom1->n_grid_1-1);i++)
	{
		e1[i]= new double[geom1->n_grid_2];
	}
	///////////////////////////////////
	//Ef//
	e2 = new double*[geom1->n_grid_1];

	for (int i=0; i<(geom1->n_grid_1);i++)
	{
		e2[i]= new double[geom1->n_grid_2];
	}
	//////////////////////////////////////
	//Ez//
	//////////////////////////////////////
	e3 = new double*[geom1->n_grid_1];

	for (int i=0; i<(geom1->n_grid_1);i++)
	{
		e3[i]= new double[geom1->n_grid_2-1];
	}
//////////////////////////////////////////
	//fi//
/////////////////////////////////////////
	fi = new double*[geom1->n_grid_1];

	for (int i=0; i<(geom1->n_grid_1);i++)
	{
		fi[i]= new double[geom1->n_grid_2];
	}
////////////////////////////////////////
	//charge_density//
///////////////////////////////////////
   t_charge_density = new double*[geom1->n_grid_1];
	   for (int i=0; i<(geom1->n_grid_1);i++)
		 {
			t_charge_density[i] = new double[geom1->n_grid_2];
		 }
/////////////////////////////////////
	   //set zero initial temp ro//
	for (int i=0; i<(geom1->n_grid_1);i++)
		for (int k=0; k<(geom1->n_grid_2);k++)
		{
			t_charge_density[i][k]=0;
		}
	
}
//Десткуктор. вивілнення памяті
E_field::~E_field()
{
	for (int i=0; i<(geom1->n_grid_1-1);i++)
		delete[]e1[i];
    delete[]e1;

	for (int i=0; i<(geom1->n_grid_1);i++)
		delete[]e2[i];
    delete[]e2;

	for (int i=0; i<(geom1->n_grid_1);i++)
		delete[]e3[i];
    delete[]e3;

	for (int i=0; i<(geom1->n_grid_1);i++)
		delete[]fi[i];
    delete[]fi;

	for (int i=0; i<(geom1->n_grid_1);i++)
	    delete[]t_charge_density[i];
    delete[]t_charge_density;
}
	//початковий розподіл E//
////////////////////////////////////////////
///////////////////////////////////////////
void E_field::set_efield_0()
{
	////Er////
	for(int i=0;i<(geom1->n_grid_1-1);i++)
		for(int k=1;k<(geom1->n_grid_2-1);k++)
		{
			e1[i][k]=0;
		}
	///Ef////
	for(int i=1;i<(geom1->n_grid_1-1);i++)
		for(int k=1;k<(geom1->n_grid_2-1);k++)
		{
			e2[i][k]=0;
		}
	///Ez////
	for(int i=1;i<(geom1->n_grid_1-1);i++)
		for(int k=0;k<(geom1->n_grid_2-1);k++)
		{
			e3[i][k]=0;
		}
}
/////////////////////////////////////////////////



////////////////////////////////////////////////

void E_field::boundary_conditions()
{
	//граничні умови Er//
////////////////////////////////
	//last elements of array [ngrid-1]!!!//
///////////////////////////////
for(int i=0;i<(geom1->n_grid_1-1);i++)
{
	e1[i][0]=0;
	e1[i][geom1->n_grid_2-1]=0;
	//e1[i][0]=limit_conditions();
}
//////////////////////////////
	
     //граничні умови Ef//
/////////////////////////////
for(int k=0;k<(geom1->n_grid_2);k++)	
{
	e2[0][k]=0;
	e2[geom1->n_grid_1-1][k]=0;
}
for(int i=0;i<(geom1->n_grid_1);i++)	
{
	e2[i][0]=0;
	e2[i][geom1->n_grid_2-1]=0;
}

////////////////////////////
	//граничні умови Ez//
/////////////////////////////////
for(int k=0;k<(geom1->n_grid_2-1);k++)	
{
	e3[0][k]=0;
	e3[geom1->n_grid_1-1][k]=0;
	//e1[0][k]=limit_conditions();
}
////////////////////////////////
}

//Electric field calculation with absorbing fields on walls//
///////////////////////////////////////////////////////////////////////////////
void E_field::calc_field(H_field* h_field1, Time* time1, current* current1, PML* pml1)
{

	double epsilon0=8.85E-12;
    int i=0;
	int k=0;
	double** j1= current1->get_j1();
	double** j2= current1->get_j2();
	double**j3= current1->get_j3();

///////Er first[i] value////
	for(k=1;k<(geom1->n_grid_2-1);k++)
		{
			 i=0;
			this->e1[i][k]=this->e1[i][k]*(2*geom1->epsilon[i][k]-geom1->sigma[i][k])/(2*geom1->epsilon[i][k]+geom1->sigma[i][k])-(j1[i][k]+(h_field1->h2[i][k]-h_field1->h2[i][k-1])/geom1->dz)*time1->delta_t/(2*geom1->epsilon[i][k]+geom1->sigma[i][k])/epsilon0;
		}
/////////////////////////////////////
	//Ez=- j on axis//
    for(k=1;k<(geom1->n_grid_2-1);k++)
		{
			 i=0;
			this->e3[i][k]=this->e3[i][k]*(2.0*geom1->epsilon[i][k]-geom1->sigma[i][k])/(2.0*geom1->epsilon[i][k]+geom1->sigma[i][k])-(j3[i][k])*time1->delta_t/(2.0*geom1->epsilon[i][k]+geom1->sigma[i][k])/epsilon0;
		}

	for( i=1;i<(geom1->n_grid_1-1);i++)
		for( k=1;k<(geom1->n_grid_2-1);k++)
		{
			this->e1[i][k]=this->e1[i][k]*(2.0*geom1->epsilon[i][k]-geom1->sigma[i][k])/(2.0*geom1->epsilon[i][k]+geom1->sigma[i][k])-(j1[i][k]+(h_field1->h2[i][k]-h_field1->h2[i][k-1])/geom1->dz)*time1->delta_t/(2.0*geom1->epsilon[i][k]+geom1->sigma[i][k])/epsilon0;
			this->e2[i][k]=this->e2[i][k]*(2.0*geom1->epsilon[i][k]-geom1->sigma[i][k])/(2.0*geom1->epsilon[i][k]+geom1->sigma[i][k])-(j2[i][k]-(h_field1->h1[i][k]-h_field1->h1[i][k-1])/geom1->dz + (h_field1->h3[i][k]-h_field1->h3[i-1][k])/geom1->dr)*time1->delta_t/(2.0*geom1->epsilon[i][k]+geom1->sigma[i][k])/epsilon0;;
			this->e3[i][k]=this->e3[i][k]*(2.0*geom1->epsilon[i][k]-geom1->sigma[i][k])/(2.0*geom1->epsilon[i][k]+geom1->sigma[i][k])-(j3[i][k]-(h_field1->h2[i][k]-h_field1->h2[i-1][k])/geom1->dr - (h_field1->h2[i][k]+h_field1->h2[i-1][k])/(2.0*geom1->dr*i))*time1->delta_t/(2.0*geom1->epsilon[i][k]+geom1->sigma[i][k])/epsilon0;
		}

    for(i=1;i<(geom1->n_grid_1-1);i++)		
		{
			 k=0;
			this->e3[i][k]=this->e3[i][k]*(2.0*geom1->epsilon[i][k]-geom1->sigma[i][k])/(2.0*geom1->epsilon[i][k]+geom1->sigma[i][k])-(j3[i][k]-(h_field1->h2[i][k]-h_field1->h2[i-1][k])/geom1->dr - (h_field1->h2[i][k]+h_field1->h2[i-1][k])/(2.0*geom1->dr*i))*time1->delta_t/(2.0*geom1->epsilon[i][k]+geom1->sigma[i][k])/epsilon0;
		}

return; 

};
///////////////////////////////////////////////////////////////////////////////////


//Electric field calculation sigma=0//
///////////////////////////////////////////////////////////////////////////////
void E_field::calc_field(H_field* h_field1, Time* time1, current* current1)
{


    int i=0;
	int k=0;
	double** j1= current1->get_j1();
	double** j2= current1->get_j2();
	double**j3= current1->get_j3();


///////Er first[i] value////
	for(k=1;k<(geom1->n_grid_2-1);k++)
		{
			 i=0;
			this->e1[i][k]=this->e1[i][k]-(j1[i][k])*time1->delta_t/epsilon0;
		}

	//Ez=- j on axis//
    for(k=1;k<(geom1->n_grid_2-1);k++)
		{
			 i=0;
			this->e3[i][k]=this->e3[i][k]-(j3[i][k])*time1->delta_t/epsilon0;
		}

	for( i=1;i<(geom1->n_grid_1-1);i++)
		for( k=1;k<(geom1->n_grid_2-1);k++)
		{
			this->e1[i][k]=this->e1[i][k]-(j1[i][k]+(h_field1->h2[i][k]-h_field1->h2[i][k-1])/geom1->dz)*time1->delta_t/epsilon0;
			this->e2[i][k]=this->e2[i][k]-(j2[i][k]-(h_field1->h1[i][k]-h_field1->h1[i][k-1])/geom1->dz + (h_field1->h3[i][k]-h_field1->h3[i-1][k])/geom1->dr)*time1->delta_t/epsilon0;;
			this->e3[i][k]=this->e3[i][k]-(j3[i][k]-(h_field1->h2[i][k]-h_field1->h2[i-1][k])/geom1->dr - (h_field1->h2[i][k]+h_field1->h2[i-1][k])/(2.0*geom1->dr*i))*time1->delta_t/epsilon0;
		}


    for(i=1;i<(geom1->n_grid_1-1);i++)		
		{
			 k=0;
			this->e3[i][k]=this->e3[i][k]-(j3[i][k]-(h_field1->h2[i][k]-h_field1->h2[i-1][k])/geom1->dr - (h_field1->h2[i][k]+h_field1->h2[i-1][k])/(2.0*geom1->dr*i))*time1->delta_t/epsilon0;
		}

return; 

};
///////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
	//set fi on r=z boundary//
/////////////////////////////////////////////////////////////////////////////
void E_field::set_fi_on_z()
	{
		for (int k=0; k<(geom1->n_grid_2);k++)
			{
				fi[geom1->n_grid_1-1][k]=0;
			}
	}

//poisson equation solving//
///////////////////////////////////////////////////////////////////////////////
void E_field::poisson_equation(Geometry* geom1, charge_density* ro1)
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
	double** ro=ro1->get_ro();

//////////////////////////////////////////////////
	//


////////////////////////////////////////////////
		for ( k=0; k<(geom1->n_grid_2);k++)
			{
  			ro[0][k]=0;
			}
		//charge_density[1][255]=1e-4;
////////////////////////////////////////////////

	//call function for cosine transform//
////////////////////////////////////////////////
////////////////////////////////////////////////

///////////////////////////////////////////////
		//copy charge_density array in to temp array//
for( i=0;i<(geom1->n_grid_1);i++)
	for( k=0;k<(geom1->n_grid_2);k++)
	{
		t_charge_density[i][k]= ro[i][k];
	}
	
	int temp=geom1->n_grid_2;
	
		for (i=0;i<geom1->n_grid_1;i++)
		{
			four1->fast_cosine_transform(t_charge_density, temp, i, false);
		}

		//sweep method//
//////////////////////////////////////////////////////////////////

	for( k=0;k<(geom1->n_grid_2);k++)
	{
			b=2.0+ pow((geom1->dr*pi*k/(geom1->dz*geom1->n_grid_2)),2);
			d=geom1->dr*geom1->dr*t_charge_density[0][k]/epsilon0;
			alpha[1]=4.0/(2.0+b);
			beta[1]=d/(2.0+b);
			for ( i=2;i<(geom1->n_grid_1-1);i++)
			  {	
			//ay-1 - by + cy+1= = -d//
				
				a = 1.0-1.0/(2.0*(i-1));
				c = 1.0+1.0/(2.0*(i-1));
				d = geom1->dr*geom1->dr*t_charge_density[i-1][k]/epsilon0;

				alpha[i] = c/(b-alpha[i-1]*a);
				beta[i]  = (d+beta[i-1]*a)/(b-alpha[i-1]*a);

			  }
			 //   i=geom1->n_grid_1;
				//a = 1.0-1.0/(2.0*(i-1));
				//c = 1.0+1.0/(2.0*(i-1));
				//d = geom1->dr*geom1->dr*charge_density[i-1][k]/epsilon0;
			   
				fi[geom1->n_grid_1-2][k]= -d-c*fi[geom1->n_grid_1-1][k]-a*beta[geom1->n_grid_1-2]/(a*alpha[geom1->n_grid_1-2]-b);

				for(i=(geom1->n_grid_1-3);i>=0;i--)
				 {
					fi[i][k]=beta[i+1]+alpha[i+1]*fi[i+1][k];
				 }
	}



//		
//for( k=0;k<(geom1->n_grid_2);k++)
//{
//	double* c = new double [geom1->n_grid_1];
//	double* d = new double [geom1->n_grid_1];
//	double* a = new double [geom1->n_grid_1];
//	double* b = new double [geom1->n_grid_1];
//	for ( i=1;i<(geom1->n_grid_1);i++)
//	{
//			a[i] = 1.0-1.0/(2.0*(i));
//			c[i] = 1.0+1.0/(2.0*(i));
//			d[i]=-geom1->dr*geom1->dr*charge_density[i][k]/epsilon0;
//			b[i]=-2.0- pow((geom1->dr*pi*k/(geom1->dz*geom1->n_grid_2)),2);
//	}
//	b[0]=-2.0- pow((geom1->dr*pi*k/(geom1->dz*geom1->n_grid_2)),2);
//    c[0] =4/(-2+b[0]) ;
//	d[0]=geom1->dr*geom1->dr*charge_density[0][k]/epsilon0;
//	d[0]=d[0]/(-2+b[0]);
//
//	for(i = 1; i < geom1->n_grid_1; i++)
//	{
//	double id = (b[i] - c[i-1] * a[i]);	/* Division by zero risk. */
//	c[i] /= id;				/* Last value calculated is redundant. */
//	d[i] = (d[i] - d[i-1] * a[i])/id;
//	}
//	fi[geom1->n_grid_1 - 1][k] = 0;
//	for(i = geom1->n_grid_1 - 2; i >= 0; i--)
//		fi[i][k] = d[i] - c[i] * fi[i + 1][k];
//
//}



	//call function for inverse cosine transform//
////////////////////////////////////////////////////////////
	for (i=0;i<geom1->n_grid_1;i++)
		{
			int temp=geom1->n_grid_2;
			four1->fast_cosine_transform(fi, temp,i, true);
		}
////////////////////////////////////////////////////////////


	//calculate electric field//
///////////////////////////////////////////////////////////
	for (i=0;i<(geom1->n_grid_1-1);i++)
		for (k=1;k<(geom1->n_grid_2);k++)
			{
				e1[i][k]=(fi[i][k]-fi[i+1][k])/geom1->dr;
			}

	for (i=0;i<(geom1->n_grid_1);i++)
		for (k=1;k<(geom1->n_grid_2-1);k++)
			{
				e3[i][k]=(fi[i][k]-fi[i][k+1])/geom1->dz;
			}
////////////////////////////////////////////////////////////

 delete [] alpha;
 delete [] beta;
}
/////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
	//function for electric field weighting//
Triple E_field::get_field(double x1, double x3)
{
	//first cell, field  weighting?????
	int i_r = (int)ceil((x1+0.5*geom1->dr)/geom1->dr)-1;
	int k_z = (int)ceil((x3+0.5*geom1->dz)/geom1->dz)-1;

	double E_r = (e1[i_r-1][k_z]+e1[i_r][k_z])/2;
	double E_z = (e3[i_r][k_z-1]+e3[i_r][k_z])/2;

	Triple components(E_r, e2[i_r][k_z], E_z);

	return components;
}
#include "Particles.h"
#include "E_field.h"
#include "H_field.h"
#include "Time.h"
#include "Triple.h"
//#include "stdafx.h"

Particles::Particles(void)
{
}
	////////конструктор//////
////////////////////////////////////////////////
Particles::Particles(char* p_name, double p_charge, double p_mass, int p_number,
					 Geometry* geom)  : geom1(geom), c_light(3.0e8), c2(9.0e16) 
{
	name = p_name;
	charge = p_charge;
	mass = p_mass;
	number = p_number;

    //allocate memory for coordinates and velocities of particles

	x1 = new double[number];
	x3 = new double[number];
	v1 = new double[number];
	v2 = new double[number];
	v3 = new double[number];
	is_alive = new bool[number];

	
};
///////////////////////////////////////////////////

//Десткуктор. вивілнення памяті
Particles::~Particles()
{

    delete [] x1;
	delete [] x3;
	delete [] v1;
	delete [] v2;
	delete [] v3;
	delete [] is_alive;
};
///////////////////////////////////////////////////////

void Particles::set_v_0()
{
	int i;
	for(i=0;i<number;i++)
	{
		v1[i] = 0.0;
		v2[i] = 1.0e5;
		v3[i] = 0.0;
	}
	
};
void Particles::set_x_0()
{
	int i;
	for(i=0;i<number;i++)
	{
		x1[i] = 0.5;
		x3[i] = 0.5;
	}
	
};
////////////////////////////////////////////////

double Particles::get_gamma(int i)
{
	return pow(1.0 - (v1[i]*v1[i] + v2[i]*v2[i] + v3[i]*v3[i])/c2,-0.5);
}



void Particles::step_v(E_field* e_fld, H_field* h_fld, Time* t)
{
	int i;
	double gamma, b1, b2, b3, e1, e2, e3, vv1, vv2, vv3;
	const double mu0 = 1e-6;
	Triple E_compon(0.0, 0.0, 0.0), B_compon(0.0, 0.0, 0.0);
	double const1 = charge*t->delta_t/2.0/mass, const2;
	for( i=0;i<number;i++)
		//if (is_alive[i])
		{
//			E_compon = e_fld->get_field(x1[i],x3[i]);
//	        B_compon = h_fld->get_field(x1[i],x3[i]);
			e1 = E_compon.first*const1;
	        e2 = E_compon.second*const1;
	        e3 = E_compon.third*const1;
	        b1 = B_compon.first*mu0*const1;
	        b2 = B_compon.second*mu0*const1;
	        b3 = B_compon.third*mu0*const1;
			//1. Multiplication by relativistic factor
			//u(n-1/2) = gamma(n-1/2)*v(n-1/)
			gamma = get_gamma(i);
			v1[i] = gamma*v1[i];
			v2[i] = gamma*v2[i];
			v3[i] = gamma*v3[i];

			//2. Half acceleration in the electric field
			//u'(n) = u(n-1/2) + q*dt/2/m*E(n)

			v1[i] = v1[i] + e1;
			v2[i] = v2[i] + e2;
			v3[i] = v3[i] + e3;

			//3. Rotation in the magnetic field
			//u" = u' + 2/(1+B'^2)[(u' + [u'xB'(n)])xB'(n)]
			//B'(n) = B(n)*q*dt/2/mass/gamma(n)
			gamma = get_gamma(i);
			b1 = b1/gamma;
	        b2 = b2/gamma;
	        b3 = b3/gamma;
			const2 = 2.0/(1.0 + b1*b1 + b2*b2 + b3*b3);
			vv1 = v1[i];
			vv2 = v2[i];
			vv3 = v3[i];
			v1[i] = vv1 + const2*((vv2 - vv1*b3 + vv3*b1)*b3 - (vv3 + vv1*b2 - vv2*b1)*b2);
			v2[i] = vv2 + const2*(-(vv1 + vv2*b3 - vv3*b2)*b3 + (vv3 + vv1*b2 - vv2*b1)*b1); 
			v3[i] = vv3 + const2*((vv1 + vv2*b3 - vv3*b2)*b2 - (vv2 - vv1*b3 + vv3*b1)*b1);

			//4. Half acceleration in the electric field
			//u(n+1/2) = u(n) + q*dt/2/m*E(n)
			v1[i] = v1[i] + e1;
			v2[i] = v2[i] + e2;
			v3[i] = v3[i] + e3;

			//5. Division by relativistic factor
			gamma = get_gamma(i);
		    v1[i] = v1[i]/gamma;
			v2[i] = v2[i]/gamma;
			v3[i] = v3[i]/gamma;

		}

}

void Particles::half_step_coord(Time* t)
{
	int i;
	for( i=0;i<number;i++)
		if (is_alive[i])
		{	
			x1[i] = x1[i] + v1[i]*t->delta_t/2.0;
            x3[i] = x3[i] + v3[i]*t->delta_t/2.0;
		}
}

////////////////////////////////////////////////////

///////////////////////////////////////////////////
	//function for charge density weighting///
void Particles::charge_weighting(charge_density* ro1)
{
	int r_cell_number=0;
	int z_cell_number=0;
	double div= 0;
	double value = 0;

	double ro_p = charge/(geom1->dr*geom1->dz);
	for(int i=0;i<number;i++)
	{
		if((x1[i]<geom1->n_grid_1)&(x3[i]<geom1->n_grid_2))
			{
				r_cell_number = (int)ceil(x1[i]/geom1->dr)-1;
				z_cell_number =  (int)ceil(x3[i]/geom1->dz)-1;
			}
		value = ((r_cell_number+1)*geom1->dr-x1[i])*((z_cell_number+1)*geom1->dz-x3[i])*ro_p;
		ro1->set_ro_weighting(r_cell_number, z_cell_number, value);

		value = (x1[i]-r_cell_number*geom1->dr)*((z_cell_number+1)*geom1->dz-x3[i])*ro_p;
		ro1->set_ro_weighting(r_cell_number+1,z_cell_number, value);

		value = (x1[i]-r_cell_number*geom1->dr)*(x3[i]-z_cell_number*geom1->dz)*ro_p;
		ro1->set_ro_weighting(r_cell_number+1,z_cell_number+1,value);

		value = ((r_cell_number+1)*geom1->dr-x1[i])*(x3[i]-z_cell_number*geom1->dz)*ro_p;
		ro1->set_ro_weighting(r_cell_number,z_cell_number+1,value);

	}
}
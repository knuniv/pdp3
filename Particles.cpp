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
	double pi = 3.14159;

	int r_i=0;  // number of particle i cell 
	int z_k=0;  // numbr of particle k cell
	
	double dr = geom1->dr;
	double dz = geom1->dz;
	double r1, r2, r3; // temp variables for calculation
	double dz1, dz2;   // temp var.: width of k and k+1 cell 
	
	double ro_v_1 =0; // charge density Q/V, V - volume of elementary cell 
	double ro_v_2=0; // charge density in i+1 cell
	
	double value =0;


	for(int i=0;i<number;i++)
	{
            // finding number of i ank cell. example: dr = 0.5; r = 0.4; i =0
		////////////////////////////
		    r_i = (int)ceil(x1[i]/dr)-1;
			z_k =  (int)ceil(x3[i]/dz)-1;
		///////////////////////////

        // in first cell ohter alg. of ro_v calc
		if(r_i>dr/2)
		{
			///////////////////////////
			 r1 =  x1[i] - 0.5*dr;
			 r2 = (r_i+0.5)*dr;
			 r3 = x1[i] + 0.5*dr;
			 ro_v_1 = charge/(pi*dz*dr*dr*2.0*i);
			 ro_v_2 = charge/(pi*dz*dr*dr*2.0*(i+1));
			 dz1 = (z_k+1)*dz-x3[i];
			 dz2 = x3[i] - z_k*dz;
		   ///////////////////////////

			// weighting in ro[i][k] cell
			value = ro_v_1*pi*dz1*(r2*r2-r1*r1); 			
			ro1->set_ro_weighting(r_i, z_k, value);
		
			// weighting in ro[i+1][k] cell
			value = ro_v_2*pi*dz1*(r3*r3-r2*r2);
			ro1->set_ro_weighting(r_i+1,z_k, value);

			// weighting in ro[i][k+1] cell
			value = ro_v_1*pi*dz2*(r2*r2-r1*r1);
			ro1->set_ro_weighting(r_i, z_k+1, value);

			// weighting in ro[i+1][k+1] cell
			value = ro_v_2*pi*dz2*(r3*r3-r2*r2);
			ro1->set_ro_weighting(r_i+1, z_k+1, value);

		}
		else 
		{
			///////////////////////////
			 r1 =  x1[i] - 0.5*dr;
			 r2 = (r_i+0.5)*dr;
			 r3 = x1[i]+0.5*dr;
			 ro_v_1 = charge/(pi*dz*dr*dr/4.0); //volume of first cell
			 ro_v_2 = charge/(pi*dz*dr*dr*2.0*(i+1));
			 dz1 = (z_k+1)*dz-x3[i];
			 dz2 = x3[i] - z_k*dz;
		   ///////////////////////////

			// weighting in ro[i][k] cell
			value = ro_v_1*pi*dz1*(r2*r2-r1*r1); 			
			ro1->set_ro_weighting(r_i, z_k, value);
		
			// weighting in ro[i+1][k] cell
			value = ro_v_2*pi*dz1*(r3*r3-r2*r2);
			ro1->set_ro_weighting(r_i+1,z_k, value);

			// weighting in ro[i][k+1] cell
			value = ro_v_1*pi*dz2*(r2*r2-r1*r1);
			ro1->set_ro_weighting(r_i, z_k+1, value);

			// weighting in ro[i+1][k+1] cell
			value = ro_v_2*pi*dz2*(r3*r3-r2*r2);
			ro1->set_ro_weighting(r_i+1, z_k+1, value);



		}
		
	}
}
#include "Particles.h"
#include "E_field.h"
#include "H_field.h"
#include "Time.h"
#include "Triple.h"
#define  pi 3.14159265358979
using namespace std;
//#include "stdafx.h"

Particles::Particles(void)
{
}
	////////constructor//////
////////////////////////////////////////////////
Particles::Particles(char* p_name, double p_charge, double p_mass, int p_number,
					 Geometry* geom)  : geom1(geom), c_light(3.0e8), c2(9.0e16) 
{
	name = p_name;
	charge = p_charge*1.6e-19;
	mass = p_mass*9.1e-31;
	number = p_number;

    //allocate memory for coordinates and velocities of particles

	x1 = new double[number];
	x3 = new double[number];
	x3 = new double[number];
	v1 = new double[number];
	v2 = new double[number];
	v3 = new double[number];
	is_alive = new bool[number];

	
};
///////////////////////////////////////////////////

//////////////////////////////////////////////////
//copy constructor//
Particles::Particles(Particles &cp_particles)
{
	name = new char[strlen(cp_particles.name)];
	strcpy(name,cp_particles.name);
	charge = cp_particles.charge;
	mass = cp_particles.mass;
	c_light = cp_particles.c_light;
	c2 = cp_particles.c2;

	x1 = new double[number];
	x3 = new double[number];
	v1 = new double[number];
	v2 = new double[number];
	v3 = new double[number];
	is_alive = new bool[number];
	for (int i=0;i<cp_particles.number;i++)
	{
		x1[i] = cp_particles.x1[i];
		x3[i] = cp_particles.x3[i];
		v1[i] = cp_particles.v1[i];
		v2[i] = cp_particles.v2[i];
		v3[i] = cp_particles.v3[i];
	}

}
/////////////////////////////////////////////////

//Десткуктор. вивілнення памяті
Particles::~Particles()
{

 //   delete [] x1;
	//delete [] x3;
	//delete [] v1;
	//delete [] v2;
	//delete [] v3;
	//delete [] is_alive;
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

double Particles::get_gamma_inv(int i)
{
	return pow((v1[i]*v1[i] + v2[i]*v2[i] + v3[i]*v3[i])/c2 + 1, 0.5);
}

void Particles::step_v(E_field *e_fld, H_field *h_fld, Time* t)
{
	int i;
	double gamma, b1, b2, b3, e1, e2, e3, vv1, vv2, vv3;
	const double mu0 = 1e-6;
	Triple E_compon(0.0, 0.0, 0.0), B_compon(0.0, 0.0, 0.0);
	double const1 = charge*t->delta_t/2.0/mass, const2;
	if (t->current_time == t->start_time) const1 = const1/2.0;
	for( i=0;i<number;i++)
		//if (is_alive[i])
		{
			E_compon = e_fld->get_field(x1[i],x3[i]);
	        B_compon = h_fld->get_field(x1[i],x3[i]);
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
			gamma = get_gamma_inv(i);
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
			gamma = get_gamma_inv(i);
		    v1[i] = v1[i]/gamma;
			v2[i] = v2[i]/gamma;
			v3[i] = v3[i]/gamma;

		}

}

void Particles::half_step_coord(Time* t)
{
	int i;
	double dr = geom1->dr;
	double dz = geom1->dz;
	double x1_wall = geom1->first_size - dr/2.0;
	double x3_wall = geom1->second_size - dz/2.0;
	double half_dr = dr/2.0;
	double half_dz = dz/2.0;
	double x1_wallX2 = x1_wall*2.0;
	double x3_wallX2 = x3_wall*2.0;
	double half_dt = t->delta_t/2.0;

	for( i=0;i<number;i++)
		if (is_alive[i])
		{	
			x1[i] = x1[i] + v1[i]*half_dt; 
            x3[i] = x3[i] + v3[i]*half_dt;

			if (x1[i] > x1_wall)
			{
				x1[i] = x1_wallX2 - x1[i];
				v1[i] = -v1[i];
			}

			if (x3[i] > x3_wall)
			{
				x3[i] = x3_wallX2 - x3[i];
				v3[i] = -v3[i];
			}

			if (x1[i] < half_dr)
			{
				x1[i] = dr - x1[i];
				v1[i] = -v1[i];
			}

			if (x3[i] < half_dz)
			{
				x3[i] = dz - x3[i];
				v3[i] = -v3[i];
			}
		}
}

////////////////////////////////////////////////////

///////////////////////////////////////////////////
	//function for charge density weighting///
void Particles::charge_weighting(charge_density* ro1)
{


	int r_i=0;  // number of particle i cell 
	int z_k=0;  // number of particle k cell
	
	double dr = geom1->dr; 
	double dz = geom1->dz;
	double r1, r2, r3; // temp variables for calculation
	double dz1, dz2;   // temp var.: width of k and k+1 cell 
	
	double ro_v =0; // charge density Q/V, V - volume of particle
	double v_1 =0; // volume of [i][k] cell
	double v_2= 0; // volume of [i+1][k] cell
	//double ro_v_2=0; // charge density in i+1 cell
	
	double value =0;
	double **temp = ro1->get_ro();
	int i;
	for(i=0;i<number;i++)
	{
            // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
		////////////////////////////
		    r_i = (int)ceil((x1[i])/dr)-1;
			z_k =  (int)ceil((x3[i])/dz)-1;
		///////////////////////////

        // in first cell other alg. of ro_v calc
		if(x1[i]>dr)
		{
			///////////////////////////
			 r1 =  x1[i] - 0.5*dr;
			 r2 = (r_i+0.5)*dr;
			 r3 = x1[i] + 0.5*dr;
			 ro_v = charge/(2.0*pi*dz*dr*x1[i]);
			 v_1 = pi*dz*dr*dr*2.0*(r_i);
			 v_2 = pi*dz*dr*dr*2.0*(r_i+1);
			 dz1 = (z_k+1)*dz-x3[i];
			 dz2 = x3[i] - z_k*dz;
		   ///////////////////////////

			// weighting in ro[i][k] cell
			value = ro_v*pi*dz1*(r2*r2-r1*r1)/v_1; 			
			ro1->set_ro_weighting(r_i, z_k, value);
		
			// weighting in ro[i+1][k] cell
			value = ro_v*pi*dz1*(r3*r3-r2*r2)/v_2;
			ro1->set_ro_weighting(r_i+1,z_k, value);

			// weighting in ro[i][k+1] cell
			value = ro_v*pi*dz2*(r2*r2-r1*r1)/v_1;
			ro1->set_ro_weighting(r_i, z_k+1, value);

			// weighting in ro[i+1][k+1] cell
			value = ro_v*pi*dz2*(r3*r3-r2*r2)/v_2;
			ro1->set_ro_weighting(r_i+1, z_k+1, value);

		}
		else 
		{
			///////////////////////////
			 r1 =  x1[i] - 0.5*dr;
			 r2 = (r_i+0.5)*dr;
			 r3 = x1[i]+0.5*dr;
			 dz1 = (z_k+1)*dz-x3[i];
			 dz2 = x3[i] - z_k*dz;
			 ro_v = charge/(2.0*pi*dz*dr*x1[i]);
			 v_1 = pi*dz*dr*dr/4.0;
			 v_2 = pi*dz*dr*dr*2.0*(r_i+1);
		   ///////////////////////////

			// weighting in ro[i][k] cell
			value = ro_v*pi*dz1*(r2*r2-r1*r1)/v_1; 			
			ro1->set_ro_weighting(r_i, z_k, value);
		
			// weighting in ro[i+1][k] cell
			value = ro_v*pi*dz1*(r3*r3-r2*r2)/v_2;
			ro1->set_ro_weighting(r_i+1,z_k, value);

			// weighting in ro[i][k+1] cell
			value = ro_v*pi*dz2*(r2*r2-r1*r1)/v_1;
			ro1->set_ro_weighting(r_i, z_k+1, value);

			// weighting in ro[i+1][k+1] cell
			value = ro_v*pi*dz2*(r3*r3-r2*r2)/v_2;
			ro1->set_ro_weighting(r_i+1, z_k+1, value);



		}
		
	}
}

////////////////////////////////////////////////////////
//function for initial (Maxwell) distribution//

void Particles::velocity_distribution(double therm_vel)
{
	
int i = 0;
int j=0;
double R =0; // number from [0;1]
double dv = 0.005; // velocity step in calculation integral
double s =0; 
double ds =0;
double* v = new double [number]; //velocity vector
double nr = sqrt(pi/2.0)* pow(therm_vel, 3);

double ss =0;

for(i=0;i<number;i++)
{
 R = (i+0.5)/number;
 double* dss= new double[100000];

 // part of numerical integral calculation//
 while (s<R*nr)
 {

	 ds = dv*j*dv*j * exp(-dv*j*dv*j / (2*therm_vel*therm_vel) ) * dv;
	 ss  = (dv*j*dv*j * exp(-(dv*j*dv*j)/(2*therm_vel*therm_vel)));
	 s = s+ds;
	 j=j+1;
 }
 v[i] = dv*j;

}
////////////////////////////////////////////////////

//random algorithm (binary)//
/////////////////////////////////////////////////

int b1 =11349;
 int  bb = b1;
  j=0;
  for (int i=0;i<=sizeof(bb)*8;i++)
  {
	   bb>>=1;
	   j=j+1;

	  if(bb==0)
	  {
		  break;
	  }

  }

double res = 0;
int pass =0;
for (int i=0;i<=j;i++)
{
     pass = b1&(1<<i);
	
	if (pass!=0)
	{
		res=res+pow(2.0,-(i+1));
	}
}		
//////////////////////////////////////////

/////////////////////////////////////////
 //random algorithm (ternary)// 
//int N =223;
//// int  bb = 0;
// int j=1;
// int aa =0;
// double res =0;
////  for (int i=4;i<=N;i++)
// // {
//	  bb = N;
//	  while(bb>=1)
//	  {
////		aa = bb % 3;
//		res = res + aa*pow(3.0, -j);
//		j=j+1;
////		bb = bb / 3;
//	  }
	 

//  }
///////////////////////////////////////////



delete []v;

}

void Particles::load_spatial_distribution(double n1, double n2)
{
	int n;
	//calculate number of electrons in a big particle
	double n_in_big = (pi*geom1->first_size*geom1->first_size*geom1->second_size/number*(n2+n1)/2.0);
	double rand_r;
	double rand_z;
	charge *= n_in_big;
	mass *= n_in_big;
	for(n = 0; n < number; n++)
	{
		rand_r = static_cast<double>(rand()) / RAND_MAX;
		rand_z = static_cast<double>(rand()) / RAND_MAX;
		x1[n] = (geom1->first_size - geom1->dr)*rand_r + geom1->dr/2.0;
		x3[n] = (geom1->second_size - geom1->dz)*sqrt(rand_z) + geom1->dz/2.0;
	}

}

void Particles::load_velocity_distribution(double v_thermal)
{
	int n;
	for (n=0; n<number; n++)
	{
		v1[n] = 0.0;
		v2[n] = 0.0;
		v3[n] = 0.0;
	}
}


void Particles::simple_j_weighting(Time* time1, current *j1, double x1_new,double x3_new, double x1_old, double x3_old, int i_n, int k_n)
{
	double dr = geom1->dr;
	double dz = geom1->dz;
	double wj = 0;
	double delta_t = time1->delta_t;

	////defining number of cell
	//int i_n = (int)ceil((x1_new)/dr)-1;
	//int k_n =(int)ceil((x3_new)/dr)-1;
	//int i_o = (int)ceil((x1_old)/dr)-1;
	//int k_o =(int)ceil((x3_old)/dr)-1;
// cheking 
	//if ((i_n !=i_o)&&(k_n!=k_o))
	//{
	//
	//}
	//else
	//{
		// distance of particle moving//
		double delta_r = x1_new - x1_old;
		double delta_z = x3_new - x3_old;

		// if i cell is not equal 0 
		if (i_n>=1)
		{
		///////////////////////////////////
		// equation y = k*x+b;//
		// finding k & b//
		double k = delta_r/delta_z;
		double b = x1_old;
	    //calculate current jz in [i,k] cell//
		wj = charge/(2*dr*dz*delta_t*2*pi*i_n*dr*dr) * (dr*delta_z - k*delta_z*delta_z/2.0 - delta_z*b + dr*dr/k * ((i_n+0.5)*(i_n+0.5)-0.25)*log((k*delta_z+b)/b)); 
		// set new weighting current value 
		j1->set_j3(i_n,k_n, wj);

        //calculate current in [i+1,k] cell//
		wj = charge/(2*dr*dz*delta_t*2*pi*(i_n+1)*dr*dr) * (k*delta_z*delta_z/2.0+delta_z*b + delta_z*dr + dr*dr/k * (0.25-(i_n+0.5)*(i_n+0.5)) * log((k*delta_z+b)/b)); 
		// set new weighting current value 
		j1->set_j3(i_n+1,k_n, wj);
		
		///////////////////////////////////
	    //calculate current jr in [i,k] cell//
		// equation y = k*x+b;//
		// finding k & b//
		 k = -delta_z/delta_r;
		 double r0 = (i_n+0.5)*dr;
		 double r1 =  x1_old;
		 b= (k_n+1.0)*dz - x3_old;

        //weighting jr in [i][k] cell
		wj = charge/(2*pi*r0*dz*dz*dr*delta_t) *(r0*k*delta_r+k/2.0 * delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b-k*(2*r0+r1))+ delta_r*(b-k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) + (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
		j1->set_j1(i_n,k_n, wj);

          b= x3_old- k_n*dz;;
         //weighting jr in [i][k+1] cell
		 wj = charge/(2*pi*r0*dz*dz*dr*delta_t) * (-r0*k*delta_r - k/2.0*delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b+k*(2*r0+r1))+ delta_r*(b+k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) - (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
		j1->set_j1(i_n, k_n+1, wj);
		}
		// if i cell is equal 0 
		else
		{
///////////////////////////////////
		// equation y = k*x+b;//
		// finding k & b//
		double k = delta_r/delta_z;
		double b = x1_old;
	    //calculate current jz in [i,k] cell//
		wj = charge/(2.0*dr*dz*delta_t*pi*dr*dr/4.0) * (dr*delta_z - k*delta_z*delta_z/2.0 - delta_z*b ); 
		// set new weighting current value 
		j1->set_j3(i_n,k_n, wj);

        //calculate current in [i+1,k] cell//
		wj = charge/(2.0*dr*dz*delta_t*2.0*pi*dr*dr) * (k*delta_z*delta_z/2.0 + delta_z*dr +delta_z*b); 
		// set new weighting current value 
		j1->set_j3(i_n+1,k_n, wj);
		
		///////////////////////////////////
	    //calculate current jr in [i,k] cell//
		// equation y = k*x+b;//
		// finding k & b//
		 k = -delta_z/delta_r;
		 double r0 = (i_n+0.5)*dr;
		 double r1 =  x1_old;
		 b= (k_n+1.0)*dz - x3_old;

        //weighting jr in [i][k] cell
		wj = charge/(2*pi*r0*dz*dz*dr*delta_t) *(r0*k*delta_r+k/2.0 * delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b-k*(2*r0+r1))+ delta_r*(b-k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) + (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
		j1->set_j1(i_n,k_n, wj);

          b= x3_old- k_n*dz;;
         //weighting jr in [i][k+1] cell
		 wj = charge/(2*pi*r0*dz*dz*dr*delta_t) * (-r0*k*delta_r - k/2.0*delta_r*(x1_old+delta_r/2.0)+0.5*delta_r*(b+k*(2*r0+r1))+ delta_r*(b+k*r1)*(4*r0*r0-dr*dr)/(8*x1_old*(x1_old+delta_r)) - (k*(r0*r0/2.0-dr*dr/8.0))*log((x1_old+delta_r)/x1_old));
		j1->set_j1(i_n, k_n+1, wj);
		}
	
	//}
}


void Particles::j_weighting(Time* time1, current *j1, Particles *old_part)
{

	double dr = geom1->dr;
	double dz = geom1->dz;

	double **J1 = j1->get_j1();
	double **J2 = j1->get_j2();
	double **J3 = j1->get_j3();
	int i;


//////////////////////////////////////////////////////////////
	for (i=0;i<number;i++)
	{
	    double x1_old=old_part->x1[i];
		double x3_old = old_part->x3[i];
		//finding number new and old cells
		int i_n = (int)ceil((x1[i])/dr)-1;
		int k_n =(int)ceil((x3[i])/dr)-1;
		int i_o = (int)ceil((x1_old)/dr)-1;
		int k_o =(int)ceil((x3_old)/dr)-1;
		if (x1_old==(i_o+1)*dr)
			i_o=i_o+1;
		if(x3_old==(k_o+1)*dz)
			k_o=k_o+1;
	    int res_cell = abs(i_n-i_o) + abs(k_n-k_o); 
		if ((x1[i]==x1_old)||(x3[i]==x3_old))
		{
			strict_motion_weighting(time1, j1,x1[i],x3[i],x1_old,x3_old);
		}
		else
		{
			switch (res_cell)
			{
			/// 1) charge in four cells
			case 0: simple_j_weighting(time1, j1, x1[i],x3[i] ,x1_old,x3_old, i_n, k_n);
			break;

           /// 2) charge in seven cells 
			case 1:
		    {
			 /// charge in seven cells (i_new != i_old)
		     if ((i_n!=i_o)&&(k_n==k_o))
			 {
				if (x1_old >(i_n+1)*dr)
					{
						double a = (x1_old - x1[i])/(x3_old - x3[i]);
						double r_boundary = (i_n+1)*dr;
						double delta_r = r_boundary - x1[i];
						double z_boundary = x3[i] + delta_r/a;

						simple_j_weighting(time1,j1, r_boundary,z_boundary ,x1_old,x3_old, i_n+1,k_n);
	   					simple_j_weighting(time1,j1, x1[i],x3[i], r_boundary,z_boundary, i_n, k_n);
					}
				else 
					{
						double a = (x1[i] - x1_old)/(x3[i] - x3_old);
						double r_boundary = (i_n)*dr;
						double delta_r = r_boundary - x1_old;
						double z_boundary = x3_old + delta_r/a;

						simple_j_weighting(time1,j1, r_boundary,z_boundary ,x1_old,x3_old, i_n-1, k_n);
						simple_j_weighting(time1,j1, x1[i],x3[i], r_boundary,z_boundary, i_n, k_n);
					}

			 }
			//  charge in seven cells (k_new != k_old)
				else if ((i_n==i_o)&&(k_n!=k_o))
					{		
						if (x3_old<k_n*dz)
							{
								double z_boundary = k_n*dz;
								double delta_z  = z_boundary - x3_old;
								double a = (x1[i] - x1_old)/(x3[i] - x3_old);
								double r_boundary = x1_old + a*delta_z;
								simple_j_weighting(time1,j1, r_boundary,z_boundary ,x1_old,x3_old, i_n, k_n-1);
								simple_j_weighting(time1,j1, x1[i],x3[i], r_boundary,z_boundary, i_n, k_n);
							}
						else 
							{
								double z_boundary = (k_n+1)*dz;
								double delta_z  = z_boundary - x3[i];
								double a = (x1_old - x1[i])/(x3_old - x3[i]);
								double r_boundary = x1[i] + a*delta_z;
								simple_j_weighting(time1,j1, r_boundary,z_boundary ,x1_old,x3_old, i_n, k_n+1);
								simple_j_weighting(time1,j1, x1[i],x3[i], r_boundary,z_boundary,i_n, k_n);
							}
					}
		    }
	     break;

	/////////////////////////////////////////
	///////// 3) charge in 10 cells /////////
	/////////////////////////////////////////
		case 2:
	{
			/// case, when particle move from [i-1] cell to [i] cell
		if (i_o<i_n)
			/////////////////////////////////////////////////////////////////
			{
			   // case, when particle move from [i-1][k-1] -> [i][k] cell
			if(k_o<k_n)
					{
						double a = (x1[i] - x1_old)/(x3[i] - x3_old);
						double r1 = i_n*dr;
					double delta_z1 = (r1 - x1_old)/a;
						double z1 = x3_old + delta_z1;
						double z2 = k_n*dz;
						double delta_r2 = (z2-x3_old)*a;
						double r2 = x1_old+ delta_r2;
						if (z1<k_n*dz)
						{
							simple_j_weighting(time1, j1, r1, z1 ,x1_old, x3_old, i_n-1, k_n-1);
							simple_j_weighting(time1, j1, r2, z2, r1, z1, i_n, k_n-1);
							simple_j_weighting(time1, j1, x1[i], x3[i], r2, z2, i_n, k_n);
						}
						else if (z1>k_n*dz)
						{
							simple_j_weighting(time1, j1, r2, z2 ,x1_old, x3_old, i_n-1, k_n-1);
							simple_j_weighting(time1, j1, r1, z1, r2, z2,i_n-1, k_n);
							simple_j_weighting(time1, j1, x1[i], x3[i], r1, z1, i_n, k_n);
						}				
					}
			// case, when particle move from [i-1][k+1] -> [i][k] cell
				else 
					{
						double a = (x1[i] - x1_old)/(x3[i] - x3_old);
						double r1 = i_n*dr;
						double delta_z1 = (r1 - x1_old)/a;
						double z1 = x3_old + delta_z1;

						double z2 = (k_n+1)*dz;
						double delta_r2 = -(x3_old-z2)*a;
						double r2 = x1_old+ delta_r2;
						if (z1>(k_n+1)*dz)
							{
								simple_j_weighting(time1, j1, r1, z1 ,x1_old, x3_old, i_n-1, k_n+1);
								simple_j_weighting(time1, j1, r2, z2, r1, z1, i_n, k_n+1);
								simple_j_weighting(time1, j1, x1[i], x3[i], r2, z2, i_n, k_n);
							}
						else if (z1<(k_n+1)*dz)
						{
								simple_j_weighting(time1, j1, r2, z2 ,x1_old, x3_old, i_n-1, k_n+1);
								simple_j_weighting(time1, j1, r1, z1, r2, z2,i_n-1, k_n);
								simple_j_weighting(time1, j1, x1[i], x3[i], r1, z1, i_n, k_n);
							}
					}
			}
	////////////////////////////////////////////////////////////////////////////
		/// case, when particle move from [i+1] cell to [i] cell
			else if (i_o>i_n)
				{
				// case, when particle move from [i+1][k-1] -> [i][k] cell
	 			if(k_o<k_n)
					{
						double a = (x1[i] - x1_old)/(x3[i] - x3_old);
						double r1 = (i_n+1)*dr;
					double delta_z1 = -(x1_old-r1)/a;
						double z1 = x3_old + delta_z1;

						double z2 = k_n*dz;
						double delta_r2 = -(z2-x3_old)*a;
					double r2 = x1_old- delta_r2;
					
						if (z1<(k_n)*dz)
						{
							simple_j_weighting(time1, j1, r1, z1 ,x1_old, x3_old, i_n+1, k_n-1);
							simple_j_weighting(time1,j1, r2, z2, r1, z1, i_n, k_n-1);
							simple_j_weighting(time1,j1, x1[i], x3[i], r2, z2, i_n, k_n);
						}
				       else if (z1>(k_n)*dz)
						{
							simple_j_weighting(time1, j1, r2, z2 ,x1_old, x3_old, i_n+1, k_n-1);
							simple_j_weighting(time1, j1, r1, z1, r2, z2,i_n+1, k_n);
							simple_j_weighting(time1, j1, x1[i], x3[i], r1, z1, i_n, k_n);
						}

					}
			    // case, when particle move from [i+1][k+1] -> [i][k] cell
			    else if (k_o>k_n)
				 {
					double a = (x1_old-x1[i])/(x3_old-x3[i]);
					double r1 = (i_n+1)*dr;
					double delta_z1 = (r1-x1[i])/a;
					double z1 = x3[i] + delta_z1;

					double z2 = (k_n+1)*dz;
					double delta_r2 = (z2-x3[i])*a;
					double r2 = x1[i] + delta_r2;
				
					if (z1>(k_n+1)*dz)
						{
							simple_j_weighting(time1, j1, r1, z1 ,x1_old,x3_old, i_n+1, k_n+1);
							simple_j_weighting(time1, j1, r2, z2, r1, z1, i_n, k_n+1);
							simple_j_weighting(time1, j1, x1[i], x3[i], r2, z2, i_n, k_n);
						}
				    else if (z1<(k_n+1)*dz)
						{
							simple_j_weighting(time1, j1, r2, z2 ,x1_old, x3_old, i_n+1, k_n+1);
							simple_j_weighting(time1, j1, r1, z1, r2, z2,i_n+1, k_n);
							simple_j_weighting(time1, j1, x1[i], x3[i], r1, z1, i_n, k_n);
						}
		     }
			 }
       	} //close case;
        break;
      } //close switch
	} //close condition
  }//close cycle

}
void Particles::azimuthal_j_weighting(Time* time1, current *this_j)
{
  
	int r_i=0;  // number of particle i cell 
	int z_k=0;  // number of particle k cell
	
	double dr = geom1->dr; 
	double dz = geom1->dz;
	double r1, r2, r3; // temp variables for calculation
	double dz1, dz2;   // temp var.: width of k and k+1 cell 
	
	double ro_v =0; // charge density Q/V, V - volume of particle
	double v_1 =0; // volume of [i][k] cell
	double v_2= 0; // volume of [i+1][k] cell
	//double ro_v_2=0; // charge density in i+1 cell
	
	double rho =0; //charge density in cell
	double current; // j_phi in cell
	double **temp = this_j->get_j2();
	int i;

			
	for(i=0;i<number;i++)
	{
            // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
		////////////////////////////
		    r_i = (int)ceil((x1[i])/dr)-1;
			z_k =  (int)ceil((x3[i])/dz)-1;
		///////////////////////////

        // in first cell other alg. of ro_v calc
		if(x1[i]>dr)
		{
			///////////////////////////
			 r1 =  x1[i] - 0.5*dr;
			 r2 = (r_i+0.5)*dr;
			 r3 = x1[i] + 0.5*dr;
			 ro_v = charge/(2.0*pi*dz*dr*x1[i]);
			 v_1 = pi*dz*dr*dr*2.0*(r_i);
			 v_2 = pi*dz*dr*dr*2.0*(r_i+1);
			 dz1 = (z_k+1)*dz-x3[i];
			 dz2 = x3[i] - z_k*dz;
		   ///////////////////////////

			// weighting in j[i][k] cell
			rho = ro_v*pi*dz1*(r2*r2-r1*r1)/v_1;
			current = rho*v2[i];
			this_j->set_j2(r_i, z_k, current);
		
			// weighting in j[i+1][k] cell
			rho = ro_v*pi*dz1*(r3*r3-r2*r2)/v_2;
			current = rho*v2[i];
			this_j->set_j2(r_i+1,z_k,  current);

			// weighting in j[i][k+1] cell
			rho = ro_v*pi*dz2*(r2*r2-r1*r1)/v_1;
			current = rho*v2[i];
			this_j->set_j2(r_i, z_k+1,  current);

			// weighting in j[i+1][k+1] cell
			rho = ro_v*pi*dz2*(r3*r3-r2*r2)/v_2;
			current = rho*v2[i];
			this_j->set_j2(r_i+1, z_k+1,  current);

		}
		else 
		{
			///////////////////////////
			 r1 =  x1[i] - 0.5*dr;
			 r2 = (r_i+0.5)*dr;
			 r3 = x1[i]+0.5*dr;
			 dz1 = (z_k+1)*dz-x3[i];
			 dz2 = x3[i] - z_k*dz;
			 ro_v = charge/(2.0*pi*dz*dr*x1[i]);
			 v_1 = pi*dz*dr*dr/4.0;
			 v_2 = pi*dz*dr*dr*2.0*(r_i+1);
		   ///////////////////////////

			// weighting in j[i][k] cell
			rho = ro_v*pi*dz1*(r2*r2-r1*r1)/v_1;
			current = rho*v2[i];
			this_j->set_j2(r_i, z_k,  current);
		
			// weighting in j[i+1][k] cell
			rho = ro_v*pi*dz1*(r3*r3-r2*r2)/v_2;
			current = rho*v2[i];
			this_j->set_j2(r_i+1,z_k,  current);

			// weighting in j[i][k+1] cell
			rho = ro_v*pi*dz2*(r2*r2-r1*r1)/v_1;
			current = rho*v2[i];
			this_j->set_j2(r_i, z_k+1,  current);

			// weighting in j[i+1][k+1] cell
			rho = ro_v*pi*dz2*(r3*r3-r2*r2)/v_2;
			current = rho*v2[i];
			this_j->set_j2(r_i+1, z_k+1,  current);



		}
		
	}

}

void Particles:: strict_motion_weighting(Time *time1, current *this_j, double x1_new,double x3_new, double x1_old, double x3_old)
{
	double dr = geom1->dr;
	double dz = geom1->dz;
	//defining number of cell
	int i_n = (int)ceil((x1_new)/dr)-1;
	int k_n =(int)ceil((x3_new)/dr)-1;
	int i_o = (int)ceil((x1_old)/dr)-1;
	int k_o =(int)ceil((x3_old)/dr)-1;

	//stirct axis motion
//////////////////////////////////////////
	if (x1_new == x1_old)
	{

			double r1=0, r2=0,r3=0;
			double delta_z = 0.0;
			double value_part = 2*pi*x1_new*dr*dz;
			double wj_lower =0;
			r1 = x1_new-0.5*dr;
			r2 = (i_n+0.5)*dr;
			r3 = x1_new+0.5*dr;
			if (i_n==0)
			{
			   wj_lower = charge/(time1->delta_t*pi*dr*dr/4.0) * pi*(r2*r2-r1*r1)/value_part;
			}
			else
			{
			 wj_lower = charge/(time1->delta_t*2*pi*i_n*dr*dr) * pi*(r2*r2-r1*r1)/value_part;
			}
			double wj_upper =  charge/(time1->delta_t*2*pi*(i_n+1)*dr*dr) *pi*(r3*r3-r2*r2)/value_part;
			double wj=0;
		    this_j->set_j1(i_n,k_n,0.0);
			this_j->set_j1(i_n,k_n+1,0.0);
			int res_k = k_n-k_o;
			switch(res_k)
			{
				case 0: 
						{
							delta_z = x3_new - x3_old;
							wj = wj_lower*delta_z;
							this_j->set_j3(i_n,k_n,wj);
							wj = wj_upper*delta_z;
							this_j->set_j3(i_n+1,k_n,wj);
						}
				break;

				case 1:
					{
						delta_z = k_n*dz - x3_old;
						wj = wj_lower*delta_z;
						this_j->set_j3(i_n,k_n-1,wj_lower);
						wj = wj_upper*delta_z;
						this_j->set_j3(i_n+1,k_n-1,wj_lower);

						delta_z = x3_new - k_n*dz;
						wj = wj_lower*delta_z;
						this_j->set_j3(i_n,k_n,wj);
						wj = wj_upper*delta_z;
						this_j->set_j3(i_n+1,k_n,wj);
					}
				break;

				case -1:
					{
						
						delta_z = (k_n+1)*dz - x3_old;
						wj = wj_lower*delta_z;
						this_j->set_j3(i_n,k_n+1,wj_lower);
						wj = wj_upper*delta_z;
						this_j->set_j3(i_n+1,k_n+1,wj_lower);

						delta_z = x3_new - (k_n+1)*dz;
						wj = wj_lower*delta_z;
						this_j->set_j3(i_n,k_n,wj);
						wj = wj_upper*delta_z;
						this_j->set_j3(i_n+1,k_n,wj);

					}
				break;
		}
   	}
///////////////////////////////////////////////////////
	
	////stirct radial motion///
//////////////////////////////////////////////////////
	else if (x3_new==x3_old)
	{
		double r0  =(i_n+0.5)*dr;
		double wj= 0;
		double delta_r=0;
		double left_delta_z = 0, right_delta_z = 0;
		double res_j = 0;
		int res_i = i_n - i_o;
		switch(res_i)
		{
		case  0:
			{
				 delta_r = x1_new - x1_old;
				 left_delta_z = (k_n+1)*dz-x3_new;
				 right_delta_z = x3_new - k_n*dz;
				 wj = charge/(pi*4.0*r0*dz*dz*dr*time1->delta_t)*(delta_r - r0*r0/(x1_old+delta_r) +r0*r0/x1_old + dr*dr/(4.0*(x1_old+delta_r)) - dr*dr/(4.0*x1_old));
				 res_j = wj*left_delta_z;
				 this_j->set_j1(i_n,k_n,res_j);
				 res_j = wj*right_delta_z;
				 this_j->set_j1(i_n,k_n+1,res_j);

			}
			break;
		case 1:
			{
				 delta_r = (i_n)*dr- x1_old;
				 left_delta_z = (k_n+1)*dz-x3_new;
				 right_delta_z = x3_new - k_n*dz;
				 r0 = (i_n-0.5)*dr;
				 wj = charge/(pi*4.0*r0*dz*dz*dr*time1->delta_t)*(delta_r - r0*r0/(x1_old+delta_r) +r0*r0/x1_old + dr*dr/(4.0*(x1_old+delta_r)) - dr*dr/(4.0*x1_old));
				 res_j = wj*left_delta_z;
				 this_j->set_j1(i_n-1,k_n,res_j);
				 res_j = wj*right_delta_z;
				 this_j->set_j1(i_n-1,k_n+1,res_j);
				
				 delta_r = x1_new - i_n*dr;
				 r0 = (i_n+0.5)*dr;
				 wj = charge/(pi*4.0*r0*dz*dz*dr*time1->delta_t)*(delta_r - r0*r0/(i_n*dr+delta_r) +r0*r0/i_n*dr + dr*dr/(4.0*(i_n*dr+delta_r)) - dr*dr/(4.0*i_n*dr));
				 res_j = wj*left_delta_z;
				 this_j->set_j1(i_n,k_n,res_j);
				 res_j = wj*right_delta_z;
				 this_j->set_j1(i_n,k_n+1,res_j);
			}
			break;
			case -1:
			{
			     delta_r = (i_n+1)*dr - x1_old ;
				 left_delta_z = (k_n+1)*dz-x3_new;
				 right_delta_z = x3_new - k_n*dz;
				 r0 = (i_n+1.5)*dr;
				 wj = charge/(pi*4.0*r0*dz*dz*dr*time1->delta_t)*(delta_r - r0*r0/(x1_old+delta_r) +r0*r0/x1_old + dr*dr/(4.0*(x1_old+delta_r)) - dr*dr/(4.0*x1_old));
				 res_j = wj*left_delta_z;
				 this_j->set_j1(i_n+1,k_n,res_j);
				 res_j = wj*right_delta_z;
				 this_j->set_j1(i_n+1,k_n+1,res_j);
				
				 delta_r = x1_new - (i_n+1)*dr;
				 r0 = (i_n+0.5)*dr;
				 wj = charge/(pi*4.0*r0*dz*dz*dr*time1->delta_t)*(delta_r - r0*r0/((i_n+1)*dr+delta_r) +r0*r0/(i_n+1)*dr + dr*dr/(4.0*((i_n+1)*dr+delta_r)) - dr*dr/(4.0*(i_n+1)*dr));
				 res_j = wj*left_delta_z;
				 this_j->set_j1(i_n,k_n,res_j);
				 res_j = wj*right_delta_z;
				 this_j->set_j1(i_n,k_n+1,res_j);
			}
		 break;
		}

	}
}
bool continuity_equation(Time *input_time, Geometry *input_geometry, current *input_J, charge_density *rho_old, charge_density *rho_new)
{
	double **rho_old_array = rho_old->get_ro() ;
    double **rho_new_array = rho_new->get_ro() ;
	double **J1 = input_J->get_j1() ;
	double **J3 = input_J->get_j3() ;
	double delta_rho = 1.0/(input_geometry->dz*4.0*3.1415*input_geometry->dr*input_geometry->dr) ;
	int i, k;
	bool ok = true;
	double res, tolerance = 1e-10 ;
	for (i=1;i<input_geometry->n_grid_1-1;i++)

		for (k=1;k<input_geometry->n_grid_2-1;k++)
		{
			res = (rho_new_array[i][k] - rho_old_array[i][k])/input_time->delta_t +(J3[i][k] - J3[i][k-1])/input_geometry->dz   +(J1[i][k] - J1[i-1][k])/input_geometry->dr + (J1[i][k] + J1[i-1][k])/(2.0*i*input_geometry->dr);
			if (res > tolerance)
			{
				ok = false;
				std::cout<<i<<" "<<k;
			    i = input_geometry->n_grid_1;
				break;
			}
		}
	return ok;
	/*return delta_rho -(J3[i][k] - J3[i][k-1])/input_geometry->dz  +
		(J1[i][k] + J1[i-1][k])/input_geometry->dr + (J1[i][k] + J1[i-1][k])/(i+0.5)/input_geometry->dr ;*/
}
#include "Load_init_param.h"


Load_init_param::Load_init_param(void)
{
}
Load_init_param::Load_init_param(char* xml_file_l):xml_file(xml_file_l)
{
}


Load_init_param::~Load_init_param(void)
{
}

void Load_init_param::read_xml()
{
	 
	int i=1.2e7;  
	// TiXmlDocument  xml_f("parameters.xml");
	TiXmlDocument  xml_f(xml_file);
	 xml_f.LoadFile();
	 if(xml_f.Error())
	 {
		 cout<< "can't read file";
			 return;
	 }

     TiXmlElement* root = xml_f.FirstChildElement ("Initial_parameters");
	 TiXmlElement* sub_root = root->FirstChildElement ("geometry");
	 char test [30];
	 TiXmlElement* ss_root = sub_root->FirstChildElement();
	const char* a =ss_root->GetText();
	// strcpy( test, sub_root->GetText());
//	 i = get_int_value(sub_root->FirstChildElement( ));
	 
	 cout << i;
	// cin  >> i;
	// TiXmlElement* sub_root = root->FirstChildElement ("geometry");
	// TiXmlElement* value = sub_root->FirstChildElement ();
	// double r_size  = get_double_value(value);
	 //double z_size  = get_double_value(value->NextSiblingElement());


}
int Load_init_param:: get_int_value(const char* r_str)
{
 int i =0;
 	 i = atoi(r_str);
	return i;
}
double Load_init_param:: get_double_value(const char* r_str)
{
 	double i = atof(r_str);
	return i;
}

 char* Load_init_param:: read_char(char* p_name)
{
	TiXmlDocument  xml_f(xml_file);
	 xml_f.LoadFile();
	 if(xml_f.Error())
	 {
		 cout<< "couldn't read file";
			 return 0;
	 }
	 int number= 0;
	 TiXmlElement* root = xml_f.FirstChildElement ("Initial_parameters"); 
    TiXmlElement* sub_root =root->FirstChildElement (p_name); 
	char* vul_arr = new char [50];
 	strcpy(vul_arr, sub_root->GetText());
	 
	   return vul_arr;
}

double* Load_init_param:: read_double_params(const char* p_name)
{
	TiXmlDocument  xml_f(xml_file);
	 xml_f.LoadFile();
	 if(xml_f.Error())
	 {
		 cout<< "couldn't read file";
			 return 0;
	 }
	 int number= 0;
	 TiXmlElement* root = xml_f.FirstChildElement ("Initial_parameters"); 
    TiXmlElement* sub_root =root->FirstChildElement (p_name); 
	TiXmlElement* read_elem =sub_root->FirstChildElement (); 
	  while(read_elem)
        {
       
			read_elem = read_elem->NextSiblingElement();
            number++;
        }
	  double* vul_arr  = new double [number-1];
	 for (int i=0;i<(number);i++)
	 read_elem =sub_root->FirstChildElement (); 

	   for(int i=0;i<(number);i++)
        {
		    vul_arr[i] = get_double_value(read_elem->GetText());
			read_elem=read_elem->NextSiblingElement();
        }	 
	   return vul_arr;
}
///////////////////////////////////////////
void Load_init_param:: read_load_particles()
{
	TiXmlDocument  xml_f(xml_file);
	 xml_f.LoadFile();
	 if(xml_f.Error())
	 {
		 cout<< "couldn't read file";
			 return ;
	 }
	 int number= 0;
	 TiXmlElement* root = xml_f.FirstChildElement ("Initial_parameters"); 
    TiXmlElement* sub_root =root->FirstChildElement ("Particle_Name"); 
	double* params= 0;
	 char* p_name= new char [50];
	 Particles*  prtls =0;
	  while(sub_root)
        {
       
			 const char* test = sub_root->GetText();
			 strcpy (p_name,sub_root->GetText());
			params = read_double_params(p_name);
			 prtls = new Particles(p_name,params,c_geom, p_list);
			 prtls->load_spatial_distribution(params[3],params[4],0,0);
			 prtls->velocity_distribution_v2(params[5]);
			 sub_root = sub_root->NextSiblingElement("Particle_Name");      
        }
	 
	   return;
}

////////////////////////////////////////////
void Load_init_param:: load_system()
{


	int i=0;
	//load PML parameters///
	///////////////////////////////////
	double * r_params   =  read_double_params("PML");
	c_pml = new PML(r_params);

	//load Geometry parameters///
	/////////////////////////////////////
    r_params   =  read_double_params("geometry");
    c_geom  = new Geometry (r_params ,c_pml);
	//delete r_params;
	//////////////////////////////////////////////////
	///////////////////////
	///creating field objects

	efield = new E_field(c_geom);
	hfield = new H_field(c_geom);
	efield->boundary_conditions();
	efield->set_homogeneous_efield(0.0, 0.0, 0);
	hfield->set_homogeneous_h(0.0, 0.0, 0.0);
	efield->set_fi_on_z();

	//load time parameters///
	/////////////////////////////////////
	 r_params   =  read_double_params("Time");
	
	 c_time = new Time(r_params);
	 	
	 //load particle parameters///
	/////////////////////////////////////
	 p_list = new particles_list(0); 
	 read_load_particles();
   //////////////////////////////////////

	 //Maxwell initial conditions///

 r_params   =  read_double_params("Boundary_Maxwell_conditions");
  Boundary_Maxwell_conditions maxwell_rad(efield);
maxwell_rad.specify_initial_field(c_geom,r_params[0],r_params[1],r_params[2]);

/////////////////////////////////////
// creating rho arrays
	charge_density rho_new(c_geom);
	charge_density rho_old(c_geom);

////////////////////////////////////
// boundary conditions///
char* a = read_char("Boundary_conditions");
if (get_int_value(a)==0)
{
 Fourier four1(0);
 Poisson_dirichlet dirih(c_geom);
	dirih.poisson_solve(efield, &rho_new);
}




}
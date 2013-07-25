#include<iostream>
#include "field.h"
#include "E_field.h"
#include "H_field.h"
#include "pdp3_time.h"
#include "Particles.h"
#include "Fourier.h"
#include "Poisson.h"
#include "Poisson_neumann.h"
#include "Poisson_dirichlet.h"
#include "particles_list.h"
#include <fstream>
#include <math.h>
#include "Boundary_Maxwell_conditions.h"
#include "input_output_class.h"
#include "Beam.h"
#include "Bunch.h"
#include "time.h"
#include "particles_struct.h"
#include "system_host.cuh"
#include "Load_init_param.h"
#define  pi = 3.1415926535897932;
using namespace std;
Particles_struct specie;	
#define BUILD_OPENCL
int main() 
{

	clock_t start, finish;
	flcuda time_elapsed;
	Load_init_param init_param("parameters.xml");
	init_param.read_xml();
	init_param.load_system();

	init_param.Run();



};
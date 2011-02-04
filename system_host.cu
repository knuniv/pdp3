//#include <cutil.h>				// cutil32.lib
#include <string.h>
#include "system_kern.cu"


extern Particles_struct specie;
	
extern "C"
{

int iDivUp (int a, int b)
{
   return (a % b != 0) ? (a / b + 1) : (a / b);
}
void computeNumBlocks (int numPnts, int maxThreads, int &numBlocks, int &numThreads)
{

    //numThreads = min( maxThreads, numPnts );
	numThreads = 256;
    numBlocks = iDivUp ( numPnts, numThreads );
}

bool InitCUDA(void)
{
   int count = 0;
   int i = 0;

   cudaGetDeviceCount(&count);
   if(count == 0) {
           fprintf(stderr, "There is no device.\n");
           return false;
   }

   for(i = 0; i < count; i++) {
           cudaDeviceProp prop;
           if(cudaGetDeviceProperties(&prop, i) == cudaSuccess) {
                   if(prop.major >= 1) {
                           break;
                   }
           }
   }
   if(i == count) {
           fprintf(stderr, "There is no device supporting CUDA.\n");
           return false;
   }
   cudaSetDevice(i);

   printf("CUDA initialized.\n");
   return true;
}


void SetupCUDA ()
{
  //computeNumBlocks ( params.NumParticles, 256, params.NumBlocks, params.NumThreads);			// particles
  //computeNumBlocks ( params.NumBoundaries, 256, params.NumBoundaryBlocks, params.NumBoundaryThreads);			// boundaries
 // CUDA_SAFE_CALL ( cudaMalloc ( (void**) &&cuda_specie, sizeof(Particles_struct)) );	
  cudaMemcpyToSymbol ( (char *)&cuda_specie, &specie, sizeof(Particles_struct), 0, cudaMemcpyHostToDevice ) ;
  cudaMalloc ( (void**) &x1, sizeof(flcuda)* specie.number) ;
  cudaMalloc ( (void**) &x3, sizeof(flcuda)* specie.number) ;
  cudaMalloc ( (void**) &v1, sizeof(flcuda)* specie.number) ;
  cudaMalloc ( (void**) &v2, sizeof(flcuda)* specie.number) ;
  cudaMalloc ( (void**) &v3, sizeof(flcuda)* specie.number) ;	

  cudaMalloc ( (void**) &e1, sizeof(flcuda)* (specie.grid_num1 - 1) * specie.grid_num3) ;
  cudaMalloc ( (void**) &e2, sizeof(flcuda)* specie.grid_num1 * specie.grid_num3 ) ;
  cudaMalloc ( (void**) &e3, sizeof(flcuda)* specie.grid_num1 * (specie.grid_num3 - 1) ) ;	

  cudaMalloc ( (void**) &h1, sizeof(flcuda)* specie.grid_num1 * (specie.grid_num3 - 1) ) ;
  cudaMalloc ( (void**) &h2, sizeof(flcuda)* (specie.grid_num1 - 1) * (specie.grid_num3 - 1) ) ;
  cudaMalloc ( (void**) &h3, sizeof(flcuda)* (specie.grid_num1 - 1) * specie.grid_num3) ;	



  cudaThreadSynchronize ();
}

void CopySpecie2Cuda (Particles_struct specie)
{
  cudaMemcpyToSymbol ( (char *)&cuda_specie, &specie, sizeof(Particles_struct), 0, cudaMemcpyHostToDevice ) ;
  cudaThreadSynchronize ();
}

//void TransferXVToCUDA (flcuda* CPU_x1, flcuda* CPU_x3, flcuda* CPU_v1, flcuda* CPU_v2, flcuda* CPU_v3);
//void TransferEHToCUDA (flcuda* CPU_e1, flcuda* CPU_e2, flcuda* CPU_e3, flcuda* CPU_h1, flcuda* CPU_h2, flcuda* CPU_h3);
void TransferXVFromCUDA (flcuda* CPU_x1, flcuda* CPU_x3, flcuda* CPU_v1, flcuda* CPU_v2, flcuda* CPU_v3);

void TransferXVToCUDA (flcuda* CPU_x1, flcuda* CPU_x3, flcuda* CPU_v1, flcuda* CPU_v2, flcuda* CPU_v3)
{
     cudaMemcpy (x1, CPU_x1, specie.number * sizeof(flcuda), cudaMemcpyHostToDevice ) ;
	 cudaMemcpy (x3, CPU_x3, specie.number * sizeof(flcuda), cudaMemcpyHostToDevice ) ;

	 cudaMemcpy (v1, CPU_v1, specie.number * sizeof(flcuda), cudaMemcpyHostToDevice ) ;
	 cudaMemcpy (v2, CPU_v2, specie.number * sizeof(flcuda), cudaMemcpyHostToDevice ) ;
	 cudaMemcpy (v3, CPU_v3, specie.number * sizeof(flcuda), cudaMemcpyHostToDevice ) ;
	
	cudaThreadSynchronize ();
}

void TransferEHToCUDA (flcuda* CPU_e1, flcuda* CPU_e2, flcuda* CPU_e3, flcuda* CPU_h1, flcuda* CPU_h2, flcuda* CPU_h3)
{
     cudaMemcpy (e1, CPU_e1, (specie.grid_num1 - 1) * specie.grid_num3 * sizeof(flcuda), cudaMemcpyHostToDevice ) ;
	 cudaMemcpy (e2, CPU_e2, specie.grid_num1 * specie.grid_num3 * sizeof(flcuda), cudaMemcpyHostToDevice ) ;
	 cudaMemcpy (e3, CPU_e3, specie.grid_num1 * (specie.grid_num3 - 1) * sizeof(flcuda), cudaMemcpyHostToDevice ) ;

	 cudaMemcpy (h1, CPU_h1, specie.grid_num1 * (specie.grid_num3 - 1) * sizeof(flcuda), cudaMemcpyHostToDevice ) ;
	 cudaMemcpy (h2, CPU_h2, (specie.grid_num1 - 1) * (specie.grid_num3 - 1) * sizeof(flcuda), cudaMemcpyHostToDevice ) ;
	 cudaMemcpy (h3, CPU_h3, (specie.grid_num1 - 1) * specie.grid_num3 * sizeof(flcuda), cudaMemcpyHostToDevice ) ;
	
	cudaThreadSynchronize ();
}

void TransferXVFromCUDA (flcuda* CPU_x1, flcuda* CPU_x3, flcuda* CPU_v1, flcuda* CPU_v2, flcuda* CPU_v3)
{
     cudaMemcpy (CPU_x1, x1, specie.number * sizeof(flcuda), cudaMemcpyDeviceToHost ) ;
	 cudaMemcpy (CPU_x3, x3, specie.number * sizeof(flcuda), cudaMemcpyDeviceToHost ) ;

	 cudaMemcpy (CPU_v1, v1, specie.number * sizeof(flcuda), cudaMemcpyDeviceToHost ) ;
	 cudaMemcpy (CPU_v2, v2, specie.number * sizeof(flcuda), cudaMemcpyDeviceToHost ) ;
	 cudaMemcpy (CPU_v3, v3, specie.number * sizeof(flcuda), cudaMemcpyDeviceToHost ) ;
	
	cudaThreadSynchronize ();
}

//void TransferToCUDA (Particle* CPU_ParticlesArray, int numPoints )
//{
//    CUDA_SAFE_CALL( cudaMemcpy (ParticlesArray, CPU_ParticlesArray, numPoints * sizeof(Particle), cudaMemcpyHostToDevice ) );
//	cudaThreadSynchronize ();
//}
//
//void TransferFromCUDA ( Particle* CPU_ParticlesArray, int numPoints )
//{
//	CUDA_SAFE_CALL( cudaMemcpy ( CPU_ParticlesArray, ParticlesArray, numPoints * sizeof(Particle), cudaMemcpyDeviceToHost ) );
//	cudaThreadSynchronize ();		
//}

void CUDA_StepV(flcuda dt)
{
	int numBlocks = 0, numThreads = 0;
	//computeNumBlocks (cuda_specie.number, 256, numBlocks, numThreads);
	numThreads = 256;
	numBlocks = iDivUp(cuda_specie.number, numThreads);
	numBlocks = 4096;
	StepV<<<numBlocks, numThreads>>> (x1, x3, v1, v2, v3, e1, e2, e3, h1, h2, h3, dt);	
	//CUT_CHECK_ERROR( "Kernel execution failed");
	cudaThreadSynchronize ();
}

//void TransferFromCUDA ( Particle* CPU_ParticlesArray, int numPoints )
//{
//	CUDA_SAFE_CALL( cudaMemcpy ( CPU_ParticlesArray, ParticlesArray, numPoints * sizeof(Particle), cudaMemcpyDeviceToHost ) );
//	cudaThreadSynchronize ();		
//}



//void CUDA_Advance (flcuda dt)
//{
//	Advance<<< params.NumBlocks, params.NumThreads>>> (ParticlesArray, params.TimeStep);	
//	CUT_CHECK_ERROR( "Kernel execution failed");
//	cudaThreadSynchronize ();
//}
//
//void CUDA_AdvanceCoordinates ()
//{
//	AdvanceCoordinates<<< params.NumBlocks, params.NumThreads>>> (ParticlesArray, params.TimeStep);	
//	CUT_CHECK_ERROR( "Kernel execution failed");
//	cudaThreadSynchronize ();
//}
//void CUDA_AdvanceVelocities ()
//{
//	AdvanceVelocities<<< params.NumBlocks, params.NumThreads>>> (ParticlesArray, params.TimeStep);	
//	CUT_CHECK_ERROR( "Kernel execution failed");
//	cudaThreadSynchronize ();
//}
//
//
//void CUDA_CopyParticlesToGLBuffer(float3* positions)
//{
//	CopyParticlesToGLBuffer<<< params.NumBlocks, params.NumThreads>>> (ParticlesArray, positions );
//}

}
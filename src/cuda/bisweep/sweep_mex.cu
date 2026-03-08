#include "mex.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Windows.h>
#include <string.h>
#include <thrust/device_vector.h>
#include <cooperative_groups.h>

#define N_EQ1	3
#define NUM_THREADS_PER_BLOCK 512
#define MAX_KNEADING_LENGTH 2001

__device__ void stepper(const double* y, double* dydt, const double* params)
{
    double a=params[0], b = params[1], c=params[2];

    dydt[0] = a*y[1] - a*y[0] + y[1]*y[2];
    dydt[1] = c*y[0] - y[1] - y[0]*y[2];
    dydt[2] = y[0]*y[1] - b*y[2];
}

__device__ void computeFixedPointDistance(double* equilibrium, const double* params){
    double a = params[0], b = params[1], c = params[2];
    double u1, u2, u4, u5; 
        u1 = (a+c)*(a+c)-4*a;
        u2 = a + c - sqrt(u1);
        u4 = sqrt(a*b*(1-0.5*u2));
        u5 = c-0.5*u2;
        equilibrium[0] = b*u5/u4;
        equilibrium[1] = u4;
        equilibrium[2] = u5;
}

__device__ double LZ76(bool * s, int n) {
    int c=1,l=1,i=0,k=1,kmax = 1,stop=0;
    while(stop ==0) {
      if (s[i+k-1] != s[l+k-1]) {
        if (k > kmax) {
          kmax=k;
        }
        i++;
  
        if (i==l) {
          c++;
          l += kmax;
          if (l+1>n)
            stop = 1;
          else {
            i=0;
            k=1;
            kmax=1;
          }
        } else {
          k=1;
        }
      } else {
        k++;
        if (l+k > n) {
          c++;
          stop =1;
        }
      }
    }
    return double(c);
}
  

__device__ double  computePeriodNormalizedKneadingSum(bool *kneadings, unsigned kneadingsLength, unsigned periodLength ){
    double kneadingSum=0, minPeriodSum=0, currPeriodSum=0;
    unsigned i=0, normalizedPeriodIndex=0;

    if(periodLength<kneadingsLength){
        for(i=0; i<periodLength; i++) {
            currPeriodSum=0;
            for(unsigned j=0; j<periodLength; j++) {
                currPeriodSum+= 2*currPeriodSum + kneadings[i+j];
            }
            if(minPeriodSum==0 || currPeriodSum < minPeriodSum) {
                minPeriodSum = currPeriodSum;
                normalizedPeriodIndex=i;
            }
        }
    }
    //filling kneading sequence with normalized period
    for(i=0; i<periodLength; i++) {
        kneadingSum = kneadingSum + kneadings[normalizedPeriodIndex+periodLength-1-i%periodLength]*pow(2,i);
    }
    return (double)(kneadingSum)/pow(2,periodLength);
}

__device__ double getNormalizedPeriodAndKneadingSum(bool* kneadings, unsigned kneadingsLength){
    //After a long transient when the periodic orbits if any have already been reached, this method computes a kneading sum that is invariant between different cyclic permutations of the period

    bool periodFound=true;
    unsigned periodLength = kneadingsLength;
    for(unsigned currPeriod=1; currPeriod < kneadingsLength/2; currPeriod++) {
        periodFound=true;
        //Check if the kneading sequence has a period with periodicity of currPeriod
        for(unsigned i=currPeriod; i < kneadingsLength-currPeriod; i+=currPeriod) {
            for ( unsigned j=0; j<currPeriod; j++) {
                if(kneadings[j] != kneadings[i+j]) {
                    periodFound=false;
                    break;
                }
            }
            if(!periodFound) {
                break;
            }
        }
        //compute kneadingSum based on period found, if any
        if(periodFound) {
            //currPeriod is the period of the kneading sequence. So this will be normalized to a sequence with the sorted period(0's followed by 1's)
            periodLength = currPeriod;
            //return computePeriodNormalizedKneadingSum( kneadings, kneadingsLength, periodLength);
            break;
        }
       
    }
    return periodLength==kneadingsLength? LZ76(kneadings, kneadingsLength) : computePeriodNormalizedKneadingSum( kneadings, kneadingsLength, periodLength);
    //return periodLength==kneadingsLength?1:0;
}

__device__ double integrator_rk4(
    double* y_current, const double* params, const double dt, 
    const unsigned N, const unsigned stride,
    const unsigned kneadingsStart, const unsigned kneadingsEnd)
{
	unsigned i, j, k, kneadingIndex=0,  kneadingArrayIndex=0;
	double dt2, dt6;
	double y1[N_EQ1], y2[N_EQ1], k1[N_EQ1], k2[N_EQ1], k3[N_EQ1], k4[N_EQ1];
	double firstDerivativeCurrent[N_EQ1],firstDerivativePrevious;
	double kneadingsWeightedSum=0;
    double fixedPointDistance[N_EQ1];
    bool kneadings[MAX_KNEADING_LENGTH];

    bool isPreviousEventOne=false, isDerviative2FirstTimeOnThisSidePositive=0, isPreviousDerivative2Positive=0, isCurrentDerivate2Positive=0;
    bool isDerviative1FirstTimeOnThisSidePositive=0, isPreviousDerivative1Positive=0, isCurrentDerivate1Positive=0;

	dt2 = dt/2.; dt6 = dt/6.;
    computeFixedPointDistance(fixedPointDistance, params);

	for(i=1; i<N; i++)
	{

        // double current_dt = dt;
        // dt2 = current_dt/2.0;
        // dt6 = current_dt/6.0;

		for(j=0; j<stride; j++)
		{
			stepper(y_current, k1, params);
			for(k=0; k<N_EQ1; k++) y1[k] = y_current[k]+k1[k]*dt2;
			stepper(y1, k2, params);
			for(k=0; k<N_EQ1; k++) y2[k] = y_current[k]+k2[k]*dt2;
			stepper(y2, k3, params);
			for(k=0; k<N_EQ1; k++) y2[k] = y_current[k]+k3[k]*dt;
			stepper(y2, k4, params);

			//Copy latest value into y_current
			for(k=0; k<N_EQ1; k++) y_current[k] += dt6*(k1[k]+2.*(k2[k]+k3[k])+k4[k]);
		}

        stepper(y_current, firstDerivativeCurrent, params);

    if((y_current[0]-fixedPointDistance[0])*(y_current[0]-fixedPointDistance[0])+(y_current[1]-fixedPointDistance[1])*(y_current[1]-fixedPointDistance[1])+(y_current[2]-fixedPointDistance[2])*(y_current[2]-fixedPointDistance[2])<1e-4){
        return -0.5;
        break;
    }else if((y_current[0]+fixedPointDistance[0])*(y_current[0]+fixedPointDistance[0])+(y_current[1]+fixedPointDistance[1])*(y_current[1]+fixedPointDistance[1])+(y_current[2]-fixedPointDistance[2])*(y_current[2]-fixedPointDistance[2])<1e-4){
        return -1;
        break;
    }else{

    if(y_current[2]<fixedPointDistance[2]){

        if(firstDerivativePrevious * firstDerivativeCurrent[2] < 0) {

            if(y_current[0] > 0){

                        if(kneadingIndex>=kneadingsStart){

                            kneadings[kneadingArrayIndex++]=1;

                        }

                        kneadingIndex++;

            } else {

                        if(kneadingIndex>=kneadingsStart){

                            kneadings[kneadingArrayIndex++]=0;

                        }

                        kneadingIndex++;

            }
        }
    
    }

        firstDerivativePrevious = firstDerivativeCurrent[2];

		if(kneadingIndex>kneadingsEnd)

			return getNormalizedPeriodAndKneadingSum(kneadings, kneadingArrayIndex);

	    }
    }
	return -0.05;
}


__global__ void sweepThreads(double* kneadingsWeightedSumSet,
    double parameter1Start, double parameter1End, unsigned parameter1Count,
	double parameter2Start, double parameter2End, unsigned parameter2Count,
	double parameter3, unsigned whichSweep,
    double x0_initial, double y0_initial, double z0_initial,
	double dt, unsigned N, unsigned stride,
	unsigned kneadingsStart, unsigned kneadingsEnd){

    int tx=blockIdx.x * blockDim.x + threadIdx.x;

	double params[3];
	double parameter1Step = (parameter1End - parameter1Start)/(parameter1Count-1);
	double parameter2Step = (parameter2End-parameter2Start)/(parameter2Count -1);
	int i,j,k;

    if(tx<parameter1Count*parameter2Count){

        i=tx/parameter1Count;
        j=tx%parameter1Count;

        if(whichSweep==0){                          //sweep of a vs b

            params[0]=parameter1Start+i*parameter1Step;
            params[1]=parameter2Start+j*parameter2Step;
            params[2]=parameter3;

        }else if(whichSweep==1){                    //sweep of a vs c

            params[0]=parameter1Start+i*parameter1Step;
            params[1]=parameter3;
            params[2]=parameter2Start+j*parameter2Step;

        }else if(whichSweep==2){                    //sweep of b vs c

            params[0]=parameter3;
            params[1]=parameter1Start+i*parameter1Step;
            params[2]=parameter2Start+j*parameter2Step;

        }

      	// double y_initial[N_EQ1]={0+0.00000001L,0,0};
        double y_initial[N_EQ1]={x0_initial,y0_initial,z0_initial};
        // double y_initial[N_EQ1]={26,0,0};
        kneadingsWeightedSumSet[i*parameter2Count+j] = integrator_rk4(y_initial, params, dt, N, stride,
                                                                kneadingsStart, kneadingsEnd);
    }
}

/**
 * error checking routine
 */
void checkAndDisplayErrors(char *label)
{
  // we need to synchronise first to catch errors due to
  // asynchronous operations that would otherwise
  // potentially go unnoticed

  cudaError_t err;

  err = cudaDeviceSynchronize();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    mexPrintf("CUDA Error: %s (at %s)\n", e, label);  
  }

  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    mexPrintf("CUDA Error: %s (at %s)\n", e, label); 
  }
}


extern "C" {

void sweep(double* kneadingsWeightedSumSet,
        double parameter1Start, double parameter1End, unsigned parameter1Count,
		double parameter2Start, double parameter2End, unsigned parameter2Count,
		double parameter3, unsigned whichSweep,
        double x0_initial, double y0_initial, double z0_initial,
		double dt, unsigned N, unsigned stride,
		unsigned kneadingsStart, unsigned kneadingsEnd){

        int totalParameterSpaceSize = parameter1Count*parameter2Count;
        double *kneadingsWeightedSumSetGpu;

        /*allocate device memory. */
        (cudaMalloc( (void**) &kneadingsWeightedSumSetGpu, totalParameterSpaceSize*sizeof(double)));
        checkAndDisplayErrors("Memory allocation");

        /*timing*/
        cudaEvent_t start_event, stop_event;
        cudaEventCreate(&start_event) ;
        cudaEventCreate(&stop_event) ;
        cudaEventRecord(start_event, 0);

        /*Dimensions. */

        int gridXDimension = totalParameterSpaceSize/NUM_THREADS_PER_BLOCK;
        if(totalParameterSpaceSize%NUM_THREADS_PER_BLOCK!=0) {
                gridXDimension += 1;
        }
        dim3 dimGrid(gridXDimension,1);
        dim3 dimBlock(NUM_THREADS_PER_BLOCK,1);

        printf(" Num of blocks per grid:       %d\n", gridXDimension);
        printf(" Num of threads per block:     %d\n", NUM_THREADS_PER_BLOCK);
        printf(" Total Num of threads running: %d\n", gridXDimension*NUM_THREADS_PER_BLOCK);
        printf(" Parameters alphaCount=%d, lCount=%d\n",parameter1Count,parameter2Count);

        /*Call kernel(global function)*/
        sweepThreads<<<dimGrid, dimBlock>>>(kneadingsWeightedSumSetGpu, parameter1Start, parameter1End, parameter1Count, parameter2Start, parameter2End, parameter2Count, parameter3, whichSweep, x0_initial, y0_initial, z0_initial, dt, N, stride, kneadingsStart, kneadingsEnd);

        cudaDeviceSynchronize();
        cudaEventRecord(stop_event, 0);
        cudaEventSynchronize(stop_event);

        float time_kernel;
        cudaEventElapsedTime(&time_kernel, start_event, stop_event);
        printf("Total time(sec) %f\n", time_kernel/1000);

        /*copy data from device memory to memory. */
        (cudaMemcpy( kneadingsWeightedSumSet, kneadingsWeightedSumSetGpu, totalParameterSpaceSize*sizeof(double), cudaMemcpyDeviceToHost));
        checkAndDisplayErrors("Error while copying kneading sum values from device to host.");

        /*Free all allocated memory. */
        cudaFree(kneadingsWeightedSumSetGpu);
    }

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    double  parameter1Start  = mxGetScalar(prhs[0]);
    double  parameter1End    = mxGetScalar(prhs[1]);
    unsigned parameter1Count = (unsigned)mxGetScalar(prhs[2]);
    
    double  parameter2Start  = mxGetScalar(prhs[3]);
    double  parameter2End    = mxGetScalar(prhs[4]);
    unsigned parameter2Count = (unsigned)mxGetScalar(prhs[5]);
    
    double  parameter3       = mxGetScalar(prhs[6]);
    unsigned whichSweep      = (unsigned)mxGetScalar(prhs[7]);

    double x0_initial        = mxGetScalar(prhs[8]);
    double y0_initial        = mxGetScalar(prhs[9]);
    double z0_initial        = mxGetScalar(prhs[10]);
    
    double  dt               = mxGetScalar(prhs[11]);
    unsigned N               = (unsigned)mxGetScalar(prhs[12]);
    unsigned stride          = (unsigned)mxGetScalar(prhs[13]);
    
    unsigned kneadingsStart  = (unsigned)mxGetScalar(prhs[14]);
    unsigned kneadingsEnd    = (unsigned)mxGetScalar(prhs[15]);

    plhs[0] = mxCreateDoubleMatrix(parameter1Count, parameter2Count, mxREAL);
    double *matlab_output = mxGetPr(plhs[0]);

    sweep(matlab_output,
        parameter1Start, parameter1End, parameter1Count,
        parameter2Start, parameter2End, parameter2Count,
        parameter3, whichSweep,
        x0_initial, y0_initial, z0_initial,
        dt, N, stride,
        kneadingsStart, kneadingsEnd);

    cudaError_t cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        mexErrMsgIdAndTxt("MexFunction:CUDA", "Kernel error: %s", cudaGetErrorString(cudaStatus));
    }

    cudaDeviceReset();
}

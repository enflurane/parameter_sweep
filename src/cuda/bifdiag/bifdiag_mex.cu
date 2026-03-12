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
#define MAX_PARAM2_COUNT 401

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

 __device__ void integrator_rk4(
    double* y_current, const double* params, const double dt, 
    const unsigned N, const unsigned stride,
    const unsigned kneadingsStart, const unsigned kneadingsEnd, double* poincare_values)
{
    unsigned i, j, k, kneadingIndex = 0, kneadingArrayIndex = 0;
    double dt2 = dt / 2., dt6 = dt / 6.;
    double y1[N_EQ1], y2[N_EQ1], k1[N_EQ1], k2[N_EQ1], k3[N_EQ1], k4[N_EQ1];
    double DerivativeCurrent, DerivativePrevious = 0;
    double fixedPointDistance[N_EQ1];
    double Previous_y_current[N_EQ1];
    double slope1;
    double slope2;

    computeFixedPointDistance(fixedPointDistance, params);

    slope2 = fixedPointDistance[1]/fixedPointDistance[0];

    for (i = 1; i < N; i++)
    {
        for (j = 0; j < stride; j++)
        {
            stepper(y_current, k1, params);
            for (k = 0; k < N_EQ1; k++) y1[k] = y_current[k] + k1[k] * dt2;
            stepper(y1, k2, params);
            for (k = 0; k < N_EQ1; k++) y2[k] = y_current[k] + k2[k] * dt2;
            stepper(y2, k3, params);
            for (k = 0; k < N_EQ1; k++) y2[k] = y_current[k] + k3[k] * dt;
            stepper(y2, k4, params);

            for (k = 0; k < N_EQ1; k++) y_current[k] += dt6 * (k1[k] + 2. * (k2[k] + k3[k]) + k4[k]);
        }

        DerivativeCurrent = fixedPointDistance[0] * y_current[1] - fixedPointDistance[1] * y_current[0];

        if (DerivativePrevious * DerivativeCurrent < 0) {

                if (kneadingIndex >= kneadingsStart) {

                    slope1 = (Previous_y_current[1] - y_current[1]) / (Previous_y_current[0] - y_current[0]);

                    poincare_values[kneadingArrayIndex++] = abs((y_current[1] - slope1 * y_current[0]) / (slope2 - slope1));

                    if(kneadingIndex >= kneadingsEnd)

                        break;

                }

            kneadingIndex++;

        }

        DerivativePrevious = DerivativeCurrent;
        for (k = 0; k < N_EQ1; k++) Previous_y_current[k] = y_current[k];

    }

}

__global__ void sweepThreads(
    double* kneadingsWeightedSumSet,
    double parameter1Start, double parameter1End, unsigned parameter1Count,
    double parameter2Start, double parameter2End, unsigned parameter2Count,
    double parameter2, double parameter3, unsigned whichSweep,
    double x0_init, double y0_init, double z0_init,
    double dt, unsigned N, unsigned stride,
    unsigned kneadingsStart, unsigned kneadingsEnd
){
    int tx = blockIdx.x * blockDim.x + threadIdx.x;

    unsigned parameter2Step = (unsigned)floor(sqrt((double)parameter1Count));
    if (parameter2Step == 0) parameter2Step = 1;

    int i = tx / parameter2Step;
    int j = tx % parameter2Step;
    int flat_index = i * parameter2Step + j;

    if (flat_index >= parameter1Count) return;

    double parameter1Step = (parameter1End - parameter1Start) / (parameter1Count - 1);

    double params[3];
    if (whichSweep == 0) {
        params[0] = parameter1Start + flat_index * parameter1Step;
        params[1] = parameter2;
        params[2] = parameter3;
    } else if (whichSweep == 1) {
        params[0] = parameter2;
        params[1] = parameter1Start + flat_index * parameter1Step;
        params[2] = parameter3;
    } else if (whichSweep == 2) {
        params[0] = parameter2;
        params[1] = parameter3;
        params[2] = parameter1Start + flat_index * parameter1Step;
    }

    double fixedPointDistance[N_EQ1];
    computeFixedPointDistance(fixedPointDistance, params);

    double x0 = fixedPointDistance[0];
    double y0 = fixedPointDistance[1] + 1;
    double z0 = fixedPointDistance[2];
    double y_initial[N_EQ1] = {x0, y0, z0};

    double kneadingvalues[MAX_PARAM2_COUNT];
    integrator_rk4(y_initial, params, dt, N, stride, kneadingsStart, kneadingsEnd, kneadingvalues);

    for (int k = 0; k < MAX_PARAM2_COUNT; k++) {
        kneadingsWeightedSumSet[flat_index * MAX_PARAM2_COUNT + k] = kneadingvalues[k];
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
		double parameter2, double parameter3, unsigned whichSweep,
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
        sweepThreads<<<dimGrid, dimBlock>>>(kneadingsWeightedSumSetGpu, parameter1Start, parameter1End, parameter1Count, parameter2Start, parameter2End, parameter2Count, parameter2 , parameter3, whichSweep, x0_initial, y0_initial, z0_initial, dt, N, stride, kneadingsStart, kneadingsEnd);

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
    
    double  parameter2       = mxGetScalar(prhs[6]); 
    double  parameter3       = mxGetScalar(prhs[7]);
    unsigned whichSweep      = (unsigned)mxGetScalar(prhs[8]);

    double x0_initial        = mxGetScalar(prhs[9]);
    double y0_initial        = mxGetScalar(prhs[10]);
    double z0_initial        = mxGetScalar(prhs[11]);
    
    double  dt               = mxGetScalar(prhs[12]);
    unsigned N               = (unsigned)mxGetScalar(prhs[13]);
    unsigned stride          = (unsigned)mxGetScalar(prhs[14]);
    
    unsigned kneadingsStart  = (unsigned)mxGetScalar(prhs[15]);
    unsigned kneadingsEnd    = (unsigned)mxGetScalar(prhs[16]);

    plhs[0] = mxCreateDoubleMatrix(parameter1Count, parameter2Count, mxREAL);
    double *matlab_output = mxGetPr(plhs[0]);

    sweep(matlab_output,
        parameter1Start, parameter1End, parameter1Count,
        parameter2Start, parameter2End, parameter2Count,
        parameter2, parameter3, whichSweep,
        x0_initial, y0_initial, z0_initial,
        dt, N, stride,
        kneadingsStart, kneadingsEnd);

    cudaError_t cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        mexErrMsgIdAndTxt("MexFunction:CUDA", "Kernel error: %s", cudaGetErrorString(cudaStatus));
    }

    cudaDeviceReset();
}
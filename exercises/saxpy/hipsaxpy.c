#include <stdio.h>
#include <stdlib.h>

#include <hip/hip_runtime.h>

__global__ void saxpy(int n, float a, float *x, float *y);

int main(void)
{
  // Host and device array pointers

  float *x, *y, *d_x, *d_y;

  float a;
  int i, n;
  
  double checksum;

  n = 1000000;
  a = 2.0;
  
  printf("HIP saxpy with n = %d\n", n);

  // Initialise x and y

  x = (float *) malloc(n*sizeof(float));
  y = (float *) malloc(n*sizeof(float));

  for (i=0; i < n; i++)
    {
      x[i] = (float) i;
      y[i] = 3.0*(((float) (n-1)) - x[i]);
    }

  // Start of HIP offload code from lectures
  
  hipMalloc(&d_x, n*sizeof(float));
  hipMalloc(&d_y, n*sizeof(float));
  hipMemcpy(d_x, x, n*sizeof(float), hipMemcpyHostToDevice);
  hipMemcpy(d_y, y, n*sizeof(float), hipMemcpyHostToDevice);

  saxpy<<<(n+255)/256, 256>>>(n, a, d_x, d_y);

  hipMemcpy(y, d_y, n*sizeof(float), hipMemcpyDeviceToHost);
  hipFree(d_x);
  hipFree(d_y);

  // End of HIP offload code from lectures

  // Check result

  checksum = 0.0;

  for (i=0; i < n; i++)
    {
      checksum += y[i];
    }

  printf("Checksum  = %lf\n", checksum);
  printf("Reference = %lf\n", ((3.0+a)/2.0)*((double) n)*((double) (n-1)));

  free(x);
  free(y);

  return 0;
}

// Start of HIP kernel code from lectures

__global__ void saxpy(int n, float a, float *x, float *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
}

// End of HIP kernel code from lectures

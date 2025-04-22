#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

int main(void)
{
  // Array pointers

  float *x, *y;

  float a;
  int i, n;

  double checksum;

  n = 1000000;
  a = 2.0;

  printf("OpenMP saxpy with n = %d\n", n);

  // Initialise x and y

  x = (float *) malloc(n*sizeof(float));
  y = (float *) malloc(n*sizeof(float));

  for (i=0; i < n; i++)
    {
      x[i] = (float) i;
      y[i] = 3.0*(((float) (n-1)) - x[i]);
    }

  // Start of OpenMP offload code from lectures

#pragma omp target teams distribute parallel for map(tofrom:y[0:n]) map(to:a,x[0:n])
  for (int i = 0; i < n; i++)
    {
      y[i] = a * x[i] + y[i];
    }  
  
  // End of OpenMP offload code from lectures

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

#include <stdio.h>
    #include <cuda.h>
 
    int main(void)
    {
 
       int nDevices;
 
       cudaGetDeviceCount(&nDevices);
       printf("Number of CUDA devices = %d\n", nDevices);
 
       return 0;
    }

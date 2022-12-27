#include <gbwt/cuda_hello.cuh>
#include <thrust/sort.h>

__global__ void hello(char *a, int *b) {
  a[threadIdx.x] += b[threadIdx.x];
  return;
}

void cuda_hello() {
  char a[N] = "Hello \0\0\0\0\0\0";
  int b[N] = {15, 10, 6, 0, -11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  char *ad;
  int *bd;
  const int csize = N * sizeof(char);
  const int isize = N * sizeof(int);

  printf("%s", a);
  printf("\n================================================\n");
  printf("%s", a);

  cudaMalloc((void **)&ad, csize);
  cudaMalloc((void **)&bd, isize);
  cudaMemcpy(ad, a, csize, cudaMemcpyHostToDevice);
  cudaMemcpy(bd, b, isize, cudaMemcpyHostToDevice);

  dim3 dimBlock(blocksize, 1);
  dim3 dimGrid(1, 1);
  hello<<<dimGrid, dimBlock>>>(ad, bd);
  cudaMemcpy(a, ad, csize, cudaMemcpyDeviceToHost);
  cudaFree(ad);
  cudaFree(bd);

  printf("%s", a);
  printf("\n================================================\n");
  printf("before sort: %s\n", a);
  thrust::stable_sort_by_key(thrust::host, a, a + 6, a,
                             thrust::greater<char>());
  printf("after sort: %s\n", a);
  return;
}

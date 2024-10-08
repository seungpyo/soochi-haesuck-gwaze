#include <stdio.h>
#include "nr.h"
#include "nrutil.h"

void print_mat(float **m, int nrl, int nrh, int ncl, int nch)
{
  for (int i = nrl; i <= nrh; ++i)
  {
    for (int j = ncl; j <= nch; ++j)
    {
      printf("%f ", m[i][j]);
    }
    printf("\n");
  }
}

int main()
{
  char fname[32];
  for (int file_idx = 1; file_idx <= 3; ++file_idx)
  {
    sprintf(fname, "lineq_dat/lineq%d.dat", file_idx);
    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
      fprintf(stderr, "Cannot open file %s\n", fname);
      return 1;
    }
    int n, m;
    fscanf(fp, "%d %d", &n, &m);
    float **a = matrix(1, n, 1, n);
    float *b = vector(1, n);
    for (int i = 1; i <= n; ++i)
    {
      for (int j = 1; j <= n; ++j)
      {
        fscanf(fp, "%f", &a[i][j]);
      }
    }
    for (int i = 1; i <= n; ++i)
    {
      fscanf(fp, "%f", &b[i]);
    }
    fclose(fp);

    printf("Matrix A:\n");
    print_mat(a, 1, n, 1, n);
    printf("Matrix B:\n");
    print_mat(b, 1, n, 1, 1);

    gaussj(a, n, b, 1);

    printf("Matrix A after gaussj:\n");
    print_mat(a, 1, n, 1, n);
    printf("Matrix B after gaussj:\n");
    print_mat(b, 1, n, 1, 1);
  }
  return 0;
}
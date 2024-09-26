#include <stdio.h>
#include "nr.h"

#define SOLVER_TYPE_F 0
#define SOLVER_TYPE_DF 1
#define MAX_BRACKETS 100
#define MAX_NAME_LEN 128

typedef float (*F)(float);
typedef void (*X2FDF)(float, float *, float *);
typedef float (*SolverF)(float (*)(float), float, float, float);
typedef float (*SolverDF)(void (*)(float, float *, float *), float, float, float);

void solve_f(
    F f,
    SolverF solver,
    int *nb,
    float roots[])
{
  float xb1[MAX_BRACKETS], xb2[MAX_BRACKETS];
  zbrak(f, 1.0, 10.0, MAX_BRACKETS, xb1, xb2, nb);
  for (int i = 1; i <= *nb; i++)
  {
    roots[i - 1] = solver(f, xb1[i], xb2[i], 1e-6);
  }
}

void solve_df(
    F f,
    X2FDF fdf,
    SolverDF solver,
    int *nb,
    float roots[])
{
  float xb1[MAX_BRACKETS], xb2[MAX_BRACKETS];
  zbrak(f, 1.0, 10.0, MAX_BRACKETS, xb1, xb2, nb);
  for (int i = 1; i <= *nb; i++)
  {
    roots[i - 1] = solver(fdf, xb1[i], xb2[i], 1e-6);
  }
}
typedef struct Args
{
  char func_name[MAX_NAME_LEN];
  char solver_name[MAX_NAME_LEN];
  F f;
  X2FDF x2fdf;
  float x1;
  float x2;
  float xacc_mul;
  SolverF solver_f;
  SolverDF solver_df;
  short solver_type;
  int *num_brackets;
  float *roots;
  unsigned int max_brackets;
} Args;

void solve(Args *args)
{
  float xacc;
  float xb1[args->max_brackets], xb2[args->max_brackets];
  zbrak(args->f, args->x1, args->x2, args->max_brackets, xb1, xb2, args->num_brackets);
  for (int i = 1; i <= *args->num_brackets; ++i)
  {
    xacc = args->xacc_mul * (xb1[i] + xb2[i]) * 0.5;
    args->roots[i - 1] = args->solver_type == SOLVER_TYPE_F
                             ? args->solver_f(args->f, xb1[i], xb2[i], xacc)
                             : args->solver_df(args->x2fdf, xb1[i], xb2[i], xacc);
  }
}

void print_result(Args *args)
{
  printf("Function: %s\n", args->func_name);
  printf("Solver: %s\n", args->solver_name);
  printf("Found %d roots, printed in (x -> f(x)) format:\n", *args->num_brackets);
  for (int i = 0; i < *args->num_brackets; ++i)
  {
    printf("\t%f -> %f\n", args->roots[i], args->f(args->roots[i]));
  }
}

void x2fdf_bessj0(float x, float *f, float *df)
{
  *f = bessj0(x);
  *df = -bessj1(x);
}

int main()
{
  int num_brackets;
  float roots[MAX_BRACKETS];
  Args args[] = {
      {
          .f = bessj0,
          .x2fdf = x2fdf_bessj0,
          .x1 = 1.0,
          .x2 = 10.0,
          .solver_f = rtbis,
          .solver_type = SOLVER_TYPE_F,
          .func_name = "bessj0",
          .solver_name = "Bisection",
      },
      {
          .f = bessj0,
          .x2fdf = x2fdf_bessj0,
          .x1 = 1.0,
          .x2 = 10.0,
          .solver_f = rtflsp,
          .solver_type = SOLVER_TYPE_F,
          .func_name = "bessj0",
          .solver_name = "Linear interpolation",
      },
      {
          .f = bessj0,
          .x2fdf = x2fdf_bessj0,
          .x1 = 1.0,
          .x2 = 10.0,
          .solver_f = rtsec,
          .solver_type = SOLVER_TYPE_F,
          .func_name = "bessj0",
          .solver_name = "Secant",
      },
      {
          .f = bessj0,
          .x2fdf = x2fdf_bessj0,
          .x1 = 1.0,
          .x2 = 10.0,
          .solver_df = rtnewt,
          .solver_type = SOLVER_TYPE_DF,
          .func_name = "bessj0",
          .solver_name = "Newton-Raphson",
      },
      {
          .f = bessj0,
          .x2fdf = x2fdf_bessj0,
          .x1 = 1.0,
          .x2 = 10.0,
          .solver_df = rtsafe,
          .solver_type = SOLVER_TYPE_DF,
          .func_name = "bessj0",
          .solver_name = "Newton-Raphson with bracketing",
      },
  };
  for (int i = 0; i < sizeof(args) / sizeof(Args); ++i)
  {
    args[i].num_brackets = &num_brackets;
    args[i].roots = roots;
    args[i].max_brackets = MAX_BRACKETS;
    args[i].xacc_mul = 1e-6;
    solve(&args[i]);
    print_result(&args[i]);
  }
  return 0;
}
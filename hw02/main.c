#include <stdio.h>
#include <string.h>
#include "nr.h"
#include "solve.h"
#include "gravity.h"

int main()
{
  int num_brackets;
  float roots[MAX_BRACKETS];
  Args args[] = {
      {
          .solver_f = rtbis,
          .solver_type = SOLVER_TYPE_F,
          .solver_name = "Bisection",
      },
      {
          .solver_f = rtflsp,
          .solver_type = SOLVER_TYPE_F,
          .solver_name = "Linear interpolation",
      },
      {
          .solver_f = rtsec,
          .solver_type = SOLVER_TYPE_F,
          .solver_name = "Secant",
      },
      {
          .solver_df = rtnewt,
          .solver_type = SOLVER_TYPE_DF,
          .solver_name = "Newton-Raphson",
      },
      {
          .solver_df = rtsafe,
          .solver_type = SOLVER_TYPE_DF,
          .solver_name = "Newton-Raphson with bracketing",
      },
      {
          .solver_f = muller_solver,
          .solver_type = SOLVER_TYPE_F,
          .solver_name = "Muller",
      },
  };
  const int num_args = sizeof(args) / sizeof(Args);
  for (int i = 0; i < num_args; ++i)
  {
    strncpy(args[i].func_name, "Bessel J0", MAX_NAME_LEN);
    args[i].num_brackets = &num_brackets;
    args[i].roots = roots;
    args[i].max_brackets = MAX_BRACKETS;
    args[i].xacc_mul = 1e-6;
  }
  for (int i = 0; i < num_args; ++i)
  {
    args[i].f = bessj0;
    args[i].x2fdf = x2fdf_bessj0;
    args[i].x1 = 1.0;
    args[i].x2 = 10.0;
    solve(&args[i]);
    print_result(&args[i]);
  }

  Args gravityArgs = {
      .func_name = "Earth-Moon Net Gravity",
      .solver_name = "Newton-Raphson with bracketing",
      .f = earth_moon_net_gravity,
      .x2fdf = earth_moon_net_gravity_fdf,
      .x1 = earth_moon_distance_in_km * 0.01,
      .x2 = earth_moon_distance_in_km * 0.99,
      .solver_df = rtsafe,
      .solver_type = SOLVER_TYPE_DF,
      .num_brackets = &num_brackets,
      .roots = roots,
      .max_brackets = 256,
      .xacc_mul = 1e-6,
  };

  solve(&gravityArgs);
  print_result(&gravityArgs);
  return 0;
}
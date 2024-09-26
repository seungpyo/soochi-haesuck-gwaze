#ifndef _SOLVE_H
#define _SOLVE_H

#define SOLVER_TYPE_F 0
#define SOLVER_TYPE_DF 1
#define MAX_BRACKETS 100
#define MAX_NAME_LEN 128

typedef float (*F)(float x);
typedef void (*X2FDF)(float x, float *fx, float *dfx);
typedef float (*SolverF)(F func, float x1, float x2, float xacc);
typedef float (*SolverDF)(X2FDF fdf, float xmin, float xmax, float xacc);

void muller(
    F f,
    float *x0,
    float *x1,
    float *x2,
    unsigned int max_iter,
    float xacc);

float muller_solver(
    float (*func)(float),
    float x1,
    float x2,
    float xacc);
#endif /* _SOLVE_H */
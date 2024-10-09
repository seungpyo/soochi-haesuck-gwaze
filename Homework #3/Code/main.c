#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nr.h"
#include "nrutil.h"

// #define SKIP_GAUSSJ

void print_mat(float **m, int row_start, int row_end, int col_start, int col_end) {
    for (int i = row_start; i <= row_end; ++i) {
        for (int j = col_start; j <= col_end; ++j) {
            if (j > col_start) {
                printf("\t");
            }
            printf("%f ", m[i][j]);
        }
        printf("\n");
    }
}

void copy_matrix(float **src, float **dst, int row_start, int row_end, int col_start, int col_end) {
    for (int i = row_start; i <= row_end; ++i) {
        for (int j = col_start; j <= col_end; ++j) {
            dst[i][j] = src[i][j];
        }
    }
}

void flatten(float **mat, float *vec, int lo, int hi) {
    for (int i = lo; i <= hi; ++i) {
        vec[i] = mat[i][1];
    }
}

void report(const char *method, float *sol, int n) {
    printf("Method: %s, Solution: [", method);
    for (int i = 1; i <= n; ++i) {
        if (i > 1) {
            printf(", ");
        }
        printf("%f", sol[i]);
    }
    printf("]\n");
}

void gauss_jordan(float **a, float **b, int n) {
    gaussj(a, n, b, 1);
    float *gauss_x = vector(1, n);
    flatten(b, gauss_x, 1, n);
    report("Gauss-Jordan", gauss_x, n);
    free_vector(gauss_x, 1, n);
}

void lu_decomposition(float **a, float **b, int n) {
    float *bvec = vector(1, n);
    int *ivec = ivector(1, n);
    float d;
    flatten(b, bvec, 1, n);
    ludcmp(a, n, ivec, &d);
    lubksb(a, n, ivec, bvec);
    report("LU Decomposition", bvec, n);
    free_ivector(ivec, 1, n);
}

void svd(float **a, float **b, int n) {
    float *w = vector(1, n);
    float **v = matrix(1, n, 1, n);
    float *x = vector(1, n);
    float *bvec = vector(1, n);
    flatten(b, bvec, 1, n);

    svdcmp(a, n, n, w, v);

    float maxSingular = 0;
    for (int i = 1; i <= n; ++i) {
        if (w[i] > maxSingular) {
            maxSingular = w[i];
        }
    }
    for (int i = 1; i <= n; ++i) {
        if (w[i] < 1e-6 * maxSingular) {
            w[i] = 0;
        }
    }

    svbksb(a, w, v, n, n, bvec, x);
    report("SVD", x, n);

    free_vector(w, 1, n);
    free_matrix(v, 1, n, 1, n);
    free_vector(x, 1, n);
    free_vector(bvec, 1, n);
}

void improved_lu_decomposition(float **a, float **b, int n) {
    int *ivec = ivector(1, n);
    float *bvec = vector(1, n);
    float *x = vector(1, n);
    float **a_mprove = matrix(1, n, 1, n);
    float d;

    copy_matrix(a, a_mprove, 1, n, 1, n);
    flatten(b, x, 1, n);
    flatten(b, bvec, 1, n);

    ludcmp(a, n, ivec, &d);
    lubksb(a, n, ivec, x);
    report("Improved LU Decomposition", x, n);

    for (int i = 1; i <= n; ++i) {
        x[i] += ((float)rand() / (float)RAND_MAX);
    }
    report("Improved LU Decomposition (noisy)", x, n);

    mprove(a_mprove, a, n, ivec, bvec, x);
    report("Improved LU Decomposition (mprove)", x, n);

    free_vector(bvec, 1, n);
    free_ivector(ivec, 1, n);
    free_vector(x, 1, n);
    free_matrix(a_mprove, 1, n, 1, n);
}

void inverse_matrix(float **a, int n) {
    int *ivec = ivector(1, n);
    float *col = vector(1, n);
    float **y = matrix(1, n, 1, n);
    float d;

    ludcmp(a, n, ivec, &d);

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            col[j] = i == j ? 1 : 0;
        }
        lubksb(a, n, ivec, col);
        for (int j = 1; j <= n; ++j) {
            y[j][i] = col[j];
        }
    }
    printf("Inverse Matrix:\n");
    print_mat(y, 1, n, 1, n);

    free_ivector(ivec, 1, n);
    free_vector(col, 1, n);
    free_matrix(y, 1, n, 1, n);
}

void determinant(float **a, int n) {
    int *ivec = ivector(1, n);
    float d;
    ludcmp(a, n, ivec, &d);
    for (int i = 1; i <= n; ++i) {
        d *= a[i][i];
    }
    printf("Determinant: %f\n", d);
    free_ivector(ivec, 1, n);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage: %s <filename>\n", argv[0]);
        return 1;
    }
    char fname[32];
    strcpy(fname, argv[1]);
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
        fprintf(stderr, "Cannot open file %s\n", fname);
        return 1;
    }
    printf("\nFile: %s\n", fname);
    int m, n;
    fscanf(fp, "%d %d", &m, &n);
    float **a = matrix(1, m, 1, n);
    float **b = matrix(1, m, 1, 1);
    float **a_bkup = matrix(1, m, 1, n);
    float **b_bkup = matrix(1, m, 1, 1);

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            fscanf(fp, "%f", &a_bkup[i][j]);
        }
    }
    for (int i = 1; i <= m; ++i) {
        fscanf(fp, "%f", &b_bkup[i][1]);
    }
    fclose(fp);

// Gauss-Jordan
#ifndef SKIP_GAUSSJ
    copy_matrix(a_bkup, a, 1, m, 1, n);
    copy_matrix(b_bkup, b, 1, m, 1, 1);
    gauss_jordan(a, b, n);
#endif  // SKIP_GAUSSJ

    // LU Decomposition
    copy_matrix(a_bkup, a, 1, m, 1, n);
    copy_matrix(b_bkup, b, 1, m, 1, 1);
    lu_decomposition(a, b, n);

    // SVD
    copy_matrix(a_bkup, a, 1, m, 1, n);
    copy_matrix(b_bkup, b, 1, m, 1, 1);
    svd(a, b, n);

    // improved LU decomposition
    copy_matrix(a_bkup, a, 1, m, 1, n);
    copy_matrix(b_bkup, b, 1, m, 1, 1);
    improved_lu_decomposition(a, b, n);

    /////////// Inverse Matrix
    copy_matrix(a_bkup, a, 1, m, 1, n);
    inverse_matrix(a, n);

    ///////////// Determinant
    copy_matrix(a_bkup, a, 1, m, 1, n);
    determinant(a, n);

    ///////////// Clean up bye bye
    free_matrix(a, 1, m, 1, n);
    free_matrix(b, 1, m, 1, 1);
    free_matrix(a_bkup, 1, m, 1, n);
    free_matrix(b_bkup, 1, m, 1, 1);

    return 0;
}
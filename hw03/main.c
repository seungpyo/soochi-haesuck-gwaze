#include <stdio.h>
#include <stdlib.h>

#include "nr.h"
#include "nrutil.h"

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

int main() {
    char fname[32];
    for (int file_idx = 1; file_idx <= 3; ++file_idx) {
        sprintf(fname, "lineq_dat/lineq%d.dat", file_idx);
        FILE *fp = fopen(fname, "r");
        if (fp == NULL) {
            fprintf(stderr, "Cannot open file %s\n", fname);
            return 1;
        }
        printf("\nFile: %s\n", fname);
        int m, n;
        fscanf(fp, "%d %d", &m, &n);
        float **a = matrix(1, m, 1, n);
        float **b_mat = matrix(1, m, 1, 1);
        float *b = vector(1, m);
        float **c = matrix(1, n, 1, 1);
        float **a_bkup = matrix(1, m, 1, n);
        float **b_bkup = matrix(1, m, 1, 1);

        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                fscanf(fp, "%f", &a[i][j]);
            }
        }
        for (int i = 1; i <= m; ++i) {
            fscanf(fp, "%f", &b_mat[i][1]);
        }
        fclose(fp);

        copy_matrix(a, a_bkup, 1, m, 1, n);
        copy_matrix(b_mat, b_bkup, 1, m, 1, 1);

        gaussj(a, n, b_mat, 1);
        printf("Method: Gauss-Jordan, Solution: [");
        for (int i = 1; i <= n; ++i) {
            if (i > 1) {
                printf(", ");
            }
            printf("%f", b_mat[i][1]);
        }
        printf("]\n");

        copy_matrix(a_bkup, a, 1, m, 1, n);
        copy_matrix(b_bkup, b_mat, 1, m, 1, 1);

        int *indices = ivector(1, n);
        flatten(b_mat, b, 1, m);
        ludcmp(a, m, indices, c);
        lubksb(a, m, indices, b);
        printf("Method: LU Decomposition, Solution: [");
        for (int i = 1; i <= n; ++i) {
            if (i > 1) {
                printf(", ");
            }
            printf("%f", b[i]);
        }
        printf("]\n");

        copy_matrix(a_bkup, a, 1, m, 1, n);
        copy_matrix(b_bkup, b_mat, 1, m, 1, 1);
        flatten(b_mat, b, 1, m);

        // SVD
        float *w = vector(1, n);
        float **v = matrix(1, n, 1, n);
        float *x = vector(1, n);
        svdcmp(a, m, n, w, v);
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
        svbksb(a, w, v, m, n, b, x);

        printf("Method: SVD, Solution: [");
        for (int i = 1; i <= n; ++i) {
            if (i > 1) {
                printf(", ");
            }
            printf("%f", x[i]);
        }
        printf("]\n");

        free_ivector(indices, 1, n);
        free_vector(b, 1, m);
        free_matrix(a, 1, m, 1, n);
        free_matrix(b_mat, 1, m, 1, 1);
        free_matrix(c, 1, n, 1, 1);
        free_matrix(a_bkup, 1, m, 1, n);
        free_matrix(b_bkup, 1, m, 1, 1);
        free_vector(w, 1, n);
        free_matrix(v, 1, n, 1, n);
        free_vector(x, 1, n);
    }
    return 0;
}
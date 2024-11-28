#include <stdio.h>
#include <stdlib.h>

#include "nr.h"
#include "nrutil.h"

#define DIM_C 3

typedef struct {
    float x;
    float y;
    float x_;
    float y_;
} DataPoint;

void *ec_malloc(size_t size) {
    void *ptr = malloc(size);
    if (ptr == NULL) {
        perror("Failed to malloc\n");
        exit(1);
    }
    return ptr;
}

void read_data(FILE *fp, int *n, DataPoint **data) {
    *n = 0;
    for (char c = fgetc(fp); c != EOF; c = fgetc(fp)) {
        if (c == '\n') {
            (*n)++;
        }
    }
    fseek(fp, 0, SEEK_SET);

    *data = (DataPoint *)ec_malloc((*n + 1) * sizeof(DataPoint));
    for (int i = 1; i <= *n; i++) {
        fscanf(fp, "%f %f %f %f", &((*data)[i].x), &((*data)[i].y), &((*data)[i].x_), &((*data)[i].y_));
    }

    fclose(fp);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <data file>\n", argv[0]);
        exit(1);
    }
    FILE *fp = fopen(argv[1], "r");
    if (fp == NULL)
        perror("File opening error\n");

    int n, i, j, k;
    DataPoint *data;
    read_data(fp, &n, &data);

    float **Jt = (float **)ec_malloc((DIM_C + 1) * sizeof(float *));
    for (i = 1; i <= DIM_C; i++) {
        Jt[i] = (float *)ec_malloc((n + 1) * sizeof(float));
    }
    for (i = 1; i <= n; i++) {
        Jt[1][i] = data[i].x;
        Jt[2][i] = data[i].y;
        Jt[3][i] = 1.0;
    }

    float **Jty_mat = (float **)ec_malloc((DIM_C + 1) * sizeof(float *));
    for (i = 1; i <= DIM_C; i++) {
        Jty_mat[i] = (float *)ec_malloc(3 * sizeof(float));
        Jty_mat[i][1] = 0.0;
        Jty_mat[i][2] = 0.0;
        for (j = 1; j <= n; j++) {
            Jty_mat[i][1] += Jt[i][j] * data[j].x_;
            Jty_mat[i][2] += Jt[i][j] * data[j].y_;
        }
    }

    float **JtJ = (float **)ec_malloc((DIM_C + 1) * sizeof(float *));
    for (i = 1; i <= DIM_C; i++) {
        JtJ[i] = (float *)ec_malloc((DIM_C + 1) * sizeof(float));
        for (j = 1; j <= DIM_C; j++) {
            JtJ[i][j] = 0.0;
            for (k = 1; k <= n; k++) {
                JtJ[i][j] += Jt[i][k] * Jt[j][k];
            }
        }
    }

    gaussj(JtJ, DIM_C, Jty_mat, 2);

    for (i = 1; i <= DIM_C; i++)
        printf("%f ", Jty_mat[i][1]);
    for (i = 1; i <= DIM_C; i++)
        printf("%f ", Jty_mat[i][2]);
    printf("\n");

    for (int i = 1; i <= DIM_C; i++) {
        free(Jt[i]);
        free(JtJ[i]);
        free(Jty_mat[i]);
    }
    free(data);
    free(Jt);
    free(JtJ);
    free(Jty_mat);

    return 0;
}
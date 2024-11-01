#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "nr.h"

typedef float (*sampler_t)(long *);

void zzamtong(int num_samples, sampler_t sampler, long *idum, float scale, float shift) {
    for (int i = 0; i < num_samples; i++) {
        printf("%f\n", sampler(idum) * scale + shift);
    }
}

void run_uniform(int num_samples, long *idum) {
    int a = -3;
    int b = 4;
    printf("Uniform(%d, %d), num_samples=%d\n", a, b, num_samples);
    zzamtong(num_samples, ran1, idum, b - a, a);
}

void run_gaussian(int num_samples, long *idum) {
    float m = 0.5;
    float s = 1.5;
    printf("Gaussian(%.2f, %.2f), num_samples=%d\n", m, s, num_samples);
    zzamtong(num_samples, gasdev, idum, s, m);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s <num_samples> <PDF: 'u' or 'g'>\n", argv[0]);
        return 1;
    }
    int num_samples = atoi(argv[1]);
    long idum = time(NULL);
    char pdf = argv[2][0];
    if (pdf == 'u') {
        run_uniform(num_samples, &idum);
    } else if (pdf == 'g') {
        run_gaussian(num_samples, &idum);
    } else {
        printf("Invalid PDF: Expected 'u' or 'g', got '%c'\n", pdf);
        return 1;
    }
    return 0;
}
#include <stdio.h>
#include "nr.h"

void get_eps(float *eps){
	*eps = 0.f;
	//Your Implementation
}

int main(){
	int ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;
	float eps, epsneg, xmin, xmax;

	machar(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp,
			&eps, &epsneg, &xmin, &xmax);
	printf("Machine Accuracy (machar): \t%0.20f\n", eps);

	get_eps(&eps);
	printf("Machine Accuracy (get_eps): \t%0.20f\n", eps);

	return 0;
}

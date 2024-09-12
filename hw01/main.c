#include <stdio.h>
#include "nr.h"

void get_eps(float *eps)
{
	float one = 1.0f;
	float prev = 1.0f, cur = 1.0f;
	do
	{
		prev = cur;
		cur /= 2.0;
	} while (one + cur != one);

	*eps = prev;
}

void get_eps_double(double *eps)
{
	double one = 1.0;
	double prev = 1.0, cur = 1.0;
	do
	{
		prev = cur;
		cur /= 2.0;
	} while (one + cur != one);

	*eps = prev;
}

int main()
{
	int ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;
	float eps, epsneg, xmin, xmax;
	double eps_double, epsneg_double, xmin_double, xmax_double;

	machar(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp,
				 &eps, &epsneg, &xmin, &xmax);
	printf("Machine Accuracy (machar): \t%0.20f\n", eps);
	get_eps(&eps);
	printf("Machine Accuracy (get_eps): \t%0.20f\n", eps);

	machar_double(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp,
								&eps_double, &epsneg_double, &xmin_double, &xmax_double);
	printf("Machine Accuracy (machar_double): \t%0.20f\n", eps_double);

	get_eps_double(&eps_double);
	printf("Machine Accuracy (get_eps_double): \t%0.20f\n", eps_double);
	return 0;
}

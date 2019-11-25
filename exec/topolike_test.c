#include "pmctools/errorlist.h"


int main(int argc, char** argv)
{
	error *myerr = NULL, **err;
	double y[3], res, logL, chi2;
	int i, n;
	double pi = 3.141592653589793;

	/* initialize */
	err = &myerr;

	y[0] = 0.0;
	y[1] = 0.0;
	y[2] = 0.0;

	n = 1000;
	for (i=0; i<n; i++) {
		y[2] = i * 2*pi/(double)n;

		use_like_(y, &res);
		chi2 = res;
		logL = -0.5 * chi2;

		printf("%g %g %g %g %g\n", y[0], y[1], y[2], chi2, logL);
	}


	return 0;
}

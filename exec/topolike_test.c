#include "exec_helper.h"

int main((int argc, char** argv)
{
	error *myerr = NULL, **err;
	double y[3], res, logL;

	/* initialize */
	err = &myerr;

	y[0] = 0.0;
	y[1] = 0.0;
	y[2] = 0.0;

	use_like_(y, &res);

	logL = -0.5 * res;

	print("%g %g %g %g\n", y[0], y[1], y[2], logL);

	return 0;
}

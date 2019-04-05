/* ============================================================ *
 * topo.c                                                       *
 * Martin Kilbinger 2019                                        *
 * ============================================================ */

#include "topo.h"

functions_wrapper_t *init_functions_topo(error **err)
{
   functions_wrapper_t *init_func;

   init_func = init_functions_wrapper(read_from_config_topo, init_topo, likeli_topo,
                                      special_topo, NULL, print_topo, err);
   forwardError(*err, __LINE__, NULL);

   return init_func;
}

void read_from_config_topo(void **state, FILE *F, error **err)
{
	topo_state *topostate;

	topostate = malloc_err(sizeof(topo_state), err);
   forwardError(*err, __LINE__,);

	*state = topostate;
}

void init_topo(common_like *like, error **err)
{

}

double likeli_topo(common_like *like, const double *params, error **err)
{
	double y[3], res;

	y[0] = params[0];
	y[1] = params[1];
	y[2] = params[2];

	use_like(y, &res);

	/* Need to return -0.5 * chi2  = -2 ln L */

	return res;
}

special_t special_topo(void *state)
{
   topo_state *topostate;
   topostate = (topo_state*)state;
   return topostate->special;
}

void print_topo(FILE *where, void *state, error **err)
{

}
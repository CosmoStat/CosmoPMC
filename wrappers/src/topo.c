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
	config_element c = {0, 0.0, ""};
   int j;

	topostate = malloc_err(sizeof(topo_state), err);
   forwardError(*err, __LINE__,);

	CONFIG_READ_S(topostate, sspecial, s, F, c, err);
   STRING2ENUM(topostate->special, topostate->sspecial, special_t, sspecial_t, j, Nspecial_t, err);


	*state = topostate;
}

void init_topo(common_like *like, error **err)
{

}

double likeli_topo(common_like *like, const double *params, error **err)
{
	double y[3], res;
	int i;

	/* Default values */
	y[0] = 0.0;
	y[1] = 0.0;
	y[2] = 0.0;

	for (i=0; i<like->npar; i++) {
      switch (like->par[i]) {
			case p_phi_topo:
				y[0] = params[i];
				break;
			case p_cos_theta_topo:
				y[1] = params[i];
				break;
			case p_alpha_topo:
				y[2] = params[i];
				break;
			default:
				break;
		}
	}

	use_like_(y, &res);

	/* Need to return -0.5 * chi2  = ln L */

	return -0.5 * res;
}

special_t special_topo(void *state)
{
   topo_state *topostate;
   topostate = (topo_state*)state;
   return topostate->special;
}

void print_topo(FILE *where, void *state, error **err)
{
   FILE *rhere;
   topo_state *topostate;

   if (where==NULL) rhere = stdin;
   else rhere = where;

   topostate = (topo_state*)state;

	fprintf(rhere, "(s)special          = (%s)%d\n", sspecial_t(topostate->special), topostate->special);

}

/* Dummy function if topolike is not linked */

#ifndef TOPO

void use_like_(double *y, double *res)
{

}

#endif

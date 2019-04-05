/* ============================================================ *
 * topo.h							*
 * Martin Kilbinger 2019					*
 * ============================================================ */

#ifndef __TOPO_H
#define __TOPO_H


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "pmctools/errorlist.h"
#include "pmctools/mvdens.h"
#include "pmctools/maths.h"

#include "param.h"
#include "wrappers.h"
#include "types.h"
#include "stdnames.h"


/* Error codes */
#define topo_base     -1300
#define topo_unknown  -1 + topo_base



typedef struct {
  mvdens *data;
  special_t special;
} topo_state;


void use_like_(double *,double *);


functions_wrapper_t *init_functions_topo(error **err);
void read_from_config_topo(void **state, FILE *F, error **err);
void init_topo(common_like *like, error **err);
double likeli_topo(common_like *like, const double *params, error **err);
special_t special_topo(void *state);
void print_topo(FILE *where, void *state, error **err);

#endif

/* ============================================================ *
 * br.h								*
 * Martin Kilbinger						*
 * 2013								*
 * ============================================================ */


#ifndef __BR_H
#define __BR_H


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "errorlist.h"

#include "bias.h"
#include "cosmo.h"
#include "lensing.h"
#include "wrappers.h"
#include "init_wrappers.h"
#include "types.h"
#include "stdnames.h"


struct bias_data_tmp; /* bias_data defined in bias.h */
typedef struct {
   char datname_gg[1024], datname_gn[1024], datname_nn[1024], covname[1024],
      model_file[1024], sspecial[128], sdata_type[128];
   struct bias_data_tmp *data;
   double invcov;              /* Anderson-Hartlap debiasing factor *
				* for inverse covariance            */
   cosmo_bias *model;
   special_t special;
   bias_data_t data_type;
} bias_state;


functions_wrapper_t *init_functions_bias(error **err);
void read_from_config_bias(void **state, FILE *F, error **err);
void init_bias(common_like *like, error **err);
void fill_parameters_bias(cosmo_bias *model, const double *params, const par_t *par, int npar, error **err);
double likeli_bias(common_like *like, const double *params, error **err);
void print_bias(FILE *where, void *state, error **err);


#endif

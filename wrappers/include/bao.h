/* ============================================================ *
 * bao.h							*
 * Martin Kilbinger 2008					*
 * ============================================================ */

#ifndef __BAO_H
#define __BAO_H


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "pmctools/errorlist.h"
#include "pmctools/mvdens.h"
#include "pmctools/maths.h"

#include "nicaea/cmb_bao.h"

#include "param.h"
#include "wrappers.h"
#include "types.h"
#include "stdnames.h"


/* Error codes */
#define bao_base     -1200
#define bao_unknown  -1 + bao_base


/* Methods: each data type can have several methods to calculate the posterior */
typedef enum {distance_A, distance_d_z, distance_D_V_ratio} method_t;
#define smethod_t(i) ( \
 i==distance_A ? "distance_A" : \
 i==distance_d_z ? "distance_d_z" : \
 i==distance_D_V_ratio ? "distance_D_V_ratio" : \
 "" )
#define Nmethod_t 3


typedef struct {
  char smethod[128], datname[1024], model_file[1024], sspecial[128];
  method_t method;
  mvdens *data;
  double *z;
  special_t special;
  cosmo *model;
} bao_state;


functions_wrapper_t *init_functions_BAO(error **err);
void read_from_config_BAO(void **state, FILE *F, error **err);
void init_BAO(common_like *like, error **err);
double likeli_BAO(common_like *like, const double *params, error **err);
special_t special_BAO(void *state);
void print_BAO(FILE *where, void *state, error **err);

#endif

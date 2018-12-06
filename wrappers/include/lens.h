/* ============================================================ *
 * lens.h							*
 * Karim Benabed, Martin Kilbinger 2008-2010			*
 * ============================================================ */

#ifndef __LENS_H
#define __LENS_H


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "errorlist.h"

#include "nicaea/cosmo.h"
#include "nicaea/lensing.h"
#include "nicaea/lensing_3rd.h"
#include "wrappers.h"
#include "init_wrappers.h"
#include "types.h"
#include "stdnames.h"

/* Error codes */
#define lens_base          -1000
#define lens_allocate      -1 + lens_base
#define lens_sigma8        -2 + lens_base
#define lens_omegabc       -3 + lens_base
#define lens_invcov        -4 + lens_base
#define lens_cov_scaling   -5 + lens_base


typedef struct {
  char datname[1024], datname2[1024], covname[1024], slensdata[128], sdecomp_eb_filter[128],
    model_file[1024], sformat[128], sspecial[128], scov_scaling[128], covname_D[1024], covname_M[1024],
    *covname_ptr[3];
  datcov *data;
  double corr_invcov;   /* Anderson-Hartlap debiasing factor for inverse covariance */
  lensdata_t lensdata;
  decomp_eb_filter_t decomp_eb_filter;
  cosebi_info_t *cosebi_info;
  lensformat_t format;
  cov_scaling_t cov_scaling;
  special_t special;
  cosmo_lens *model;
  cosmo_3rd *model_3rd;
  order_t order;
  double a1, a2;            /* Coefficients for quadratic angular bin weighting */
  double sigma_8;
} lens_state;


functions_wrapper_t *init_functions_Lensing(error **err);
void read_from_config_Lensing(void **state, FILE *F, error **err);
void init_Lensing(common_like *like, error** err);
void fill_parameters_lens(cosmo_lens *lensmodel, cosmo_3rd *lensmodel_3rd, lens_state *lensstate, const double *params, const par_t *par, int npar, error **err);
double likeli_Lensing(common_like *like, const double *params, error **err);
special_t special_Lensing(void *state);
void print_Lensing(FILE *where, void *state, error **err);

#endif

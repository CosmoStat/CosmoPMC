/* ============================================================ *
 * halo.h							*
 * Martin Kilbinger 2009				        *
 * ============================================================ */

#ifndef __HALO_H
#define __HALO_H


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "pmctools/errorlist.h"
#include "pmctools/mvdens.h"
#include "pmctools/maths.h"

#include "wrappers.h"
#include "param.h"

#include "nicaea/cmb_bao.h"
#include "nicaea/hod.h"


typedef struct {
  char datname[1024], covname[1024], shalodata[128], model_file[1024], shalomode[128], sngal_fit_type[128],
    sspecial[128], sintconst_type[128], intconst_file[128];
  wt_t *wt;
  halodata_t halodata;
  halomode_t halomode;
  cosmo_hm *model;
  double delta, intconst, ngal, ngalerr, corr_invcov;
  ngal_fit_t ngal_fit_type;
  intconst_t intconst_type;
  special_t special;
} halo_state;


functions_wrapper_t *init_functions_GalCorr(error **err);
void read_from_config_GalCorr(void **state, FILE *F, error **err);
void init_GalCorr(common_like *like, error **err);
void fill_parameters_hm(cosmo_hm *halomodel, const double *params, const par_t *par, int npar, error **err);
double likeli_GalCorr(common_like *like, const double *params, error **err);
void print_GalCorr(FILE *where, void *state, error **err);
void write_to_config_GalCorr(FILE *OUT, const void *what, error **err);


#endif

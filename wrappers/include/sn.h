/* ============================================================ *
 * sn.h								*
 * ============================================================ */

#ifndef __SN_H
#define __SN_H

#include <sys/times.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "math.h"
#include "errorlist.h"
#include "nicaea/sn1a.h"
#include "param.h"
#include "wrappers.h"
#include "types.h"
#include "stdnames.h"


/* Error codes */
#define sn_base     -800
#define sn_allocate -1 + sn_base
#define sn_dead     -2 + sn_base
#define sn_choke    -3 + sn_base
#define sn_timeout  -4 + sn_base
#define sn_readerr  -5 + sn_base
#define sn_unknown  -6 + sn_base

typedef struct {

  char datname[1024], zAV_name[1024], datname_beta_d[1024], model_file[1024], sspecial[128],
    sdatformat[128];
  char schi2mode[128];
  chi2mode_t chi2mode;
  sndatformat_t datformat;
  double Theta2_denom[3];          /* For chi2_Theta2_denom_fixed */
  interTable *zAV;		   /* For chi2_dust */
  mvdens *data_beta_d;
  special_t special;
  cosmo_SN *model;
  int add_logdetCov;

  SnSample *sample;

} Sn_state;


functions_wrapper_t *init_functions_SNIa(error **err);
void read_from_config_SNIa(void **state, FILE *F, error **err);
void init_SNIa(common_like *like, error **err);
double likeli_SNIa(common_like *like,  const double *params, error **err);
special_t special_SN(void *state);
void print_SNIa(FILE *where, void *state, error **err);


#endif 


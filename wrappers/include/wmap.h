/* Full WMAP{3,5,7} likelihood (e.g. Dunkley et al. 2008), module 'CMB' */

#ifndef __WMAP_H
#define __WMAP_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "timexec.h"
#include "errorlist.h"
#include "param.h"
#include "math.h"
#include "nicaea/cmb_bao.h"

#ifdef DOWMAP
#include "C_wrappers.h"
#endif

/* TRUE_FALSE(0) = 'F', TRUE_FALSE(1) = 'T' */
#define TRUE_FALSE(i) (i==0 ? 'F' : 'T')


/* Indices of CAMB input parameters */
typedef enum {ci_omegac=0, ci_omegab=1, ci_omegaK=2, ci_H0=3, ci_ns=4, ci_Delta2R=5,
              ci_alphas=6, ci_w0=7, ci_tau=8, ci_Neffnu0=9, ci_Neffnumass=10,
              ci_omeganumass=11, ci_nt=12, ci_r=13, ci_kpivs=14, ci_kpivt=15} ci_t;


typedef struct {
  char scamb_path[2048], data_path[2048], Cl_SZ_file[2048], model_file[1024];
  cosmo *model;
  forked_state* camb;
  forked_state* lkl;
  double *cls;
  double *parout;
  int want_sigma_8, want_out_Cl, want_out_Pk, do_lensing, is_tau;
  double* template;
  int n_template, lmax, accurate;
  double sigma_8;
} wmap_state;


functions_wrapper_t *init_functions_CMB(error **err);
void read_from_config_CMB(void **state, FILE *F, error **err);
double *read_Cl_SZ_file(const char *Cl_SZ_file, double lmax, error **err);
void init_CMB(common_like *like, error** err);
double likeli_CMB(common_like *like, const double* params, error** err);
void deduced_CMB(const common_like *like, double *res);
void print_CMB(FILE *where, void *state, error **err);
void check_scamb(void *state, error **err);
void wmap_reparametrize(common_like *self, const double *parin, double * parout, error** err);

/* For CMB using the distance priors (Komatsu et al. 2008), module CMBDistPrior */
typedef struct {
  char datname[1024], model_file[1024], sspecial[128];
  mvdens *data;
  cosmo *model;
  special_t special;
} cmbDP_state;


functions_wrapper_t *init_functions_CMBDistPrior(error **err);
void read_from_config_CMBDistPrior(void **state, FILE *F, error **err);
void init_CMBDistPrior(common_like *like, error **err);
double likeli_CMBDistPrior(common_like *like, const double *params, error **err);
special_t special_CMBDistPrior(void *state);
void print_CMBDistPrior(FILE *where, void *state, error **err);

#ifndef DOWMAP
void* lkl_init(const char *, error **);
double compute_lkl(void* ptr, double* tt, double* te, double* ee, double* bb, 
		   int *pnlmax);
void f90_compute_lkl_(double* cl_tt, double* cl_te, double *cl_ee, double *cl_bb,
		      int *nlmax, double *like_tot);
#endif

/* Maximum wave mode l, 2000 is necessary for 0.3% precision for C_l's */
#define TT_MAX 2000

#define WMAP_NOSTATUS 0
#define WMAP_OK       1
#define WMAP_TIMEOUT -2
#define WMAP_CHOKE   -3
#define WMAP_DEAD    -4

/* errors */
#define wmap_base -2000
#define wmap_allocate  -1 + wmap_base
#define wmap_dead      -2 + wmap_base
#define wmap_choke     -3 + wmap_base
#define wmap_timeout   -4 + wmap_base
#define wmap_unknown   -5 + wmap_base
#define wmap_de_prior  -6 + wmap_base
#define wmap_ell       -7 + wmap_base
#define wmap_lmax      -8 + wmap_base

#endif

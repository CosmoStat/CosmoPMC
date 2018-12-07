/* ============================================================ *
 * exec_helper.h						*
 * Martin Kilbinger, Karim Benabed				*
 * ============================================================ */

// put here some specific functions needed by executables
// everything that is specific to the parameters 

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "pmclib/pmc.h"
#include "pmctools/errorlist.h"
#include "pmctools/maths.h"
#include "config.h"
#include "stdnames.h"
#include "types.h"
#include "par.h"
#include "init_wrappers.h"


typedef enum {avg_mean, avg_median} avg_t;
#define savg_t(i) ( \
  i==avg_mean   ? "mean"   : \
  i==avg_median ? "median" : \
  "" \
)
#define Navg_t 2


#define exec_helper_base    -11000
#define err_dirname             -1 + exec_helper_base 
#define err_median		-2 + exec_helper_base

struct par1dw {
  double par, w;
};


gsl_rng *init_random(int seed, FILE *FLOG);

void print_mean_sigma(const char *pname, const double *pmean, const double *pvar, const par_t *par, int ntot,
		      avg_t avgtype, error **err);
int par1dw_cmp(const void *a, const void *b);
double median_from_psim(double *X, double *w, short *flg, int nsamples, int ndim, int a, error **err);
void sigma_from_psim(double *X, double *weights, short *flg, int nsamples, int ndim,
                     int a, double mean, double *sigma, const double *confidence, error **err);
double *mean_sigma_from_psim(pmc_simu *psim, const char *pname, const par_t *par, avg_t avgtype, error **err);
void write_mean(double *pmean, const char *name, int ndata, void **data_extra, error **err);
void write_model_raw(double *model, int Nmodel, const char *name, error **err);
void confidence_level(double *X, double *weights, short *flg, int nsamples,
                      int ndim, int a, double val, int dir, const double *conf_xxx, double *cl, error **err);

void covariance_from_sample_all(const double *X, const double *weights, int ndim, long nsamples,
                                const double *X_plus_ded, int ndim_plus_ded, const char *dir, error **err);
void covariance_from_sample(const double *X, const double *weight, int ndim, long nsample, 
                            const char *name, const char *cinvname, error **err);
void write_header_cosmo_pmc(FILE *OUT, const par_t *par, int npar, int n_ded);
void print_step_cosmo_pmc(FILE* where, double weight, double loglkl, size_t npar, size_t n_ded, const double *params,
		   const double *params_ded);
void out_pmc_simu_cosmo_pmc(const char *name, const pmc_simu *psim, const par_t *par, double norm, error **err);

void make_and_test_dir(const char *dirname, error **err);

/* ============================================================ *
 * param.h.							*
 * Martin Kilbinger 2008					*
 * Contains stuff which is used in mkmc, pmc and some wrapper   *
 * files.							*
 * ============================================================ */

#ifndef __PARAM_H
#define __PARAM_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "config.h"
#include "nhist.h"
#include "pmctools/errorlist.h"
#include "pmctools/io.h"
#include "pmctools/maths.h"
#include "pmctools/mvdens.h"
#include "pmclib/pmc.h"
#include "pmclib/mcmc.h"

#include "types.h"
#include "wrappers.h"
#include "all_wrappers.h"
#include "stdnames.h"
#include "par.h"

/* The following should go into a separate file in wrappers/ */


/* Error codes */
#define mk_base           -500
#define mk_npar           -1 + mk_base
#define mk_negative       -2 + mk_base
#define mk_n_ded          -3 + mk_base
#define mk_omegabc        -4 + mk_base
#define mk_h100           -5 + mk_base
#define mk_w1de           -6 + mk_base
#define mk_sigma8         -7 + mk_base
#define mk_null           -8 + mk_base
#define mk_empty          -9 + mk_base
#define mk_boxmax        -10 + mk_base
#define mk_chdir         -11 + mk_base
#define mk_cwdbuf        -12 + mk_base
#define mk_scamb         -13 + mk_base
#define mk_undef         -14 + mk_base
#define mk_range         -15 + mk_base
#define mk_special_prior -16 + mk_base
#define mk_data          -17 + mk_base


/* The following is defined in pmclib. It is needed
 * here as long pmc is still used internally.
 */
#ifndef tls_cosmo_par
   #define tls_cosmo_par -17 + tls_base
#endif


/* Time-out in milliseconds */
#define DEFAULT_TIMEOUT  100000


/* MCMC starting points */
typedef enum {start_ran, start_fid, start_min, start_nul, start_max} start_t;
#define sstart_t(i) ( \
 i==start_ran ? "ran" : \
 i==start_fid ? "fid" : \
 i==start_min ? "min" : \
 i==start_nul ? "nul" : \
 i==start_max ? "max" : \
 "" )
#define Nstart_t 5


/* Initial proposal, MCMC and PMC */
typedef enum {mcmcini_fisher_inv, mcmcini_fisher, mcmcini_fisher_diag, mcmcini_previous,
	      mcmcini_diag, pmcini_fisher_rshift, pmcini_fisher_eigen, pmcini_file,
              pmcini_random_pos, pmcini_none} initial_t;
#define sinitial_t(i) ( \
 i==mcmcini_fisher_inv   ? "Fisher_inv" : \
 i==mcmcini_fisher       ? "Fisher" : \
 i==mcmcini_fisher_diag  ? "Fisher_diag" : \
 i==mcmcini_previous     ? "previous" : \
 i==mcmcini_diag         ? "diag" : \
 i==pmcini_fisher_rshift ? "fisher_rshift" : \
 i==pmcini_fisher_eigen  ? "fisher_eigen" : \
 i==pmcini_file          ? "file" : \
 i==pmcini_random_pos    ? "random_pos" : \
 i==pmcini_none          ? "none" : \
 "" )
#define Ninitial_t 10

/* What to do with dead proposal components */
typedef enum {bury, revive} dead_comp_t;
#define sdead_comp_t(i) ( \
  i==bury   ? "bury"    : \
  i==revive ? "revive"  : \
  "")
#define Ndead_comp_t 2

typedef enum {tempering_none, tempering_linear, tempering_log} tempering_t;
#define stempering_t(i) ( \
  i==tempering_none   ? "none" : \
  i==tempering_linear ? "linear" : \
  i==tempering_log    ? "log" : \
  "")
#define Ntempering_t 3

/* Basic configuration file information */
typedef struct {
  int ndata, npar, n_ded, nbinhist;
  char sdata[50], sprior[500];
  mvdens *prior;
  int nprior, *indprior;
  data_t *data;
  double *min, *max, logpr_default;
  common_like **data_extra;
  par_t *par;
  initial_t initial;
  char sinitial[50], **spar;
  double version;
  int myid, nproc;
} config_base;

/* Config info for MCMC */
typedef struct {
  config_base base;
  int nchain, ncov, ndecorr;
  double *fid;
  start_t start;
  char sstart[512];
  char *dirname;
  double fudge, fburnin, boxdiv;
} config_mcmc;

/* Config info for PMC */
typedef struct {
  config_base base;
  int ndata, nsamples, ncomp, niter, nclipw;
  char sprop[64], prop_ini_name[512], sdead_comp[64];
  double fshift, fvar, fsfinal, fmin, fmax;
  dead_comp_t dead_comp;
  int df;
  char stempering[64];
  tempering_t tempering;
  double t_min;
} config_pmc;

/* Config info for max_post and Fisher matrix */
typedef struct {
  config_base base;
  start_t start;
  double *fid, tolerance;
  char sstart[512];
} config_max;

typedef config_max config_fish;


/* To use the CONFIG_READ... macros the entries have to be in a structure */
typedef struct {
  int ncomp;
  char name[1024], deduced[1024];
  double THETA_MIN, THETA_MAX;
} config_dummy;

int checkPrintErr_and_continue(error **err, FILE *F1, FILE *F2, const double *params, int npar);
mvdens *assign_prior(char sprior[], int npar, int no_init, error **err);
void read_config_base(config_base *config, FILE *F, int no_init, error **err);
void read_config_mcmc_file(config_mcmc *config, const char *cname, error **err);
void read_config_pmc_file(config_pmc *, const char *, mix_mvdens **, const gsl_rng *, int, error **);
void read_config_max_file(config_max *config, const char *cname, error **err);
void read_config_fish_file(config_max *config, const char *cname, error **err);
void comp_shift_mean(const config_pmc *config, double *mean, const double *pos0, const gsl_rng *rng,
		     double fshift, error **err);
void read_initial_proposal_from_config(FILE *F, config_pmc *config, mix_mvdens **proposal_return,
				       const gsl_rng *rng, int no_init, error **err);
void read_extra_for_data(FILE *F, data_t data, common_like *like, int no_init, error **err);
void out_extra(FILE *OUT, data_t data, void *extra, error **err);
void out_config_base(FILE *OUT, config_base config, error **err);
void out_config_mcmc(FILE *OUT, config_mcmc config, error **err);
void out_config_pmc(FILE *OUT, config_pmc config, error **err);
void out_config_max(FILE *OUT, config_max config, error **err);
void out_config_fish(FILE *OUT, config_fish config, error **err);

void write_config_base_file(FILE *OUT, config_base config, error **err);
void write_config_pmc_file(FILE *OUT, config_pmc config, error **err);

parabox *parabox_from_config(int npar, const double *min, const double *max, error **err);

double get_default_log_prior(const double *min, const double *max, int npar, error **err);
double likelihood_log_pdf_single(void *extra, const double *x, error **err);
double posterior_log_pdf_common(config_base *config, const double *x,  error **err);
double posterior_log_pdf_common_void(void *config, const double *x, error **err);
double prior_log_pdf_special(special_t special, const par_t *par, const double *min,
			     const double *max, int npar, error **err);
void retrieve_ded(const void *extra, double *x_ded, error **err);
int retrieve_ded_single(const void *data_extra, double *x_ded, data_t data);


void out_histogram(const nd_histogram *hist, FILE *F);
void out_histogram_data(const nd_histogram *hist, FILE *F, error **err);

void free_config(config_mcmc config, error **err);

void create_1d_2d_hist(double *accstates, long naccepted, int nbinhist, int npar, const double *min,
		       const double *max, double *weights, FILE *FLOG, const char *path, error **err);
void create_1d_hist(double *accstates, long naccepted, int nbinhist, int npar, const double *min,
		    const double *max, double *weights, const char *path, error **err);
void create_2d_hist(double *accstates, long naccepted, int nbinhist, int npar, const double *min,
		    const double *max, double *weights, const char *path, error **err);

void create_sigma_8_chain(const config_mcmc *config, FILE *FLOG, error **err);
void close_sigma_8_chain(const config_mcmc *config, FILE *FLOG);

void check_npar(int npar, int nmax, error **err);
double *merge_param_param_ded(double *param, const double *param_ded, long nsample, int npar, int n_ded, error **err);

void do_nothing();


common_like *extra_init_common(data_t i, const par_t *par, int npar, error **err);

special_t get_special(data_t data, common_like *extra, error **err);

functions_wrapper_t *init_functions_Mvdens(error **err);
void read_mvdens_component(mvdens *g, FILE *F, error **err);
void read_from_config_Mvdens(void **state, FILE *f, error **err);
void init_dummy(common_like *like, error **err);
double likeli_Mvdens(common_like *like, const double *params, error **err);
functions_wrapper_t *init_functions_MixMvdens(error **err);
void read_from_config_MixMvdens(void **state, FILE *F, error **err);
double likeli_MixMvdens(common_like *like, const double *params, error **err);

void set_base_parameters(cosmo *cosmo, double Omegam, double Omegab, double Omegade, double Omeganumass,
			 double Omegac, double OmegaK, double omegam, double omegab, double omegade,
			 double omeganumass, double omegac, double omegaK, double h100,
			 int iOmegade, int iOmegaK, int iomegade, int iomegaK, error **err);
void set_base_Omegam(cosmo *cosmo, double Omegam, double Omegab, double Omegade, double Omeganumass,
		     double Omegac, double OmegaK, int iOmegaK);
void set_base_Omegab(cosmo *cosmo, double Omegam, double Omegab, double Omegac);
void set_base_Omegade(cosmo *cosmo, double Omegam, double Omegab, double Omegade, double Omeganumass,
		      double Omegac, double OmegaK);
void set_base_Omeganumass(cosmo *cosmo, double Omeganumass);
void reset_base_parameters(double *Omegam, double *Omegab, double *Omegade, double *Omeganumass,
			   double *Omegac, double *OmegaK,
			   double *omegam, double *omegab, double *omegade, double *omeganumass,
			   double *omegac, double *omegaK, double *h100,
			   int *iOmegade, int *iOmegaK, int *iomegade, int *iomegaK);


#endif
 


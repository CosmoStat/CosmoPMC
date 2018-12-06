#ifndef __COSMO_MCMC_H
#define __COSMO_MCMC_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <sys/times.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "io.h"
#include "maths.h"
#include "errorlist.h"
#include "pmclib/parabox.h"
#include "config.h"
#include "mvdens.h"

//#include "cosmo.h"
//#include "sn1a.h"

#include "param.h"
#include "pmclib/mcmc.h"
#include "pmclib/tools.h"
#include "stdnames.h"

#include "wrappers.h"
#include "all_wrappers.h"

#include "mkmax.h"
#include "exec_helper.h"


/* Output covariance at each update (if applicable) */
#define OUT_COV 1


int add_prior(mvdens *A, double factor, FILE *FLOG, error **err);
int de_singularize(double *A, int n, double diag_val, FILE *F, error **err);
void print_matrix(const double *a, int n, FILE *F);
mvdens *initial_proposal(config_mcmc config, FILE *FLOG, FILE *CHFROM, double correct, error **err);
int mcmc_one(double *pstate, double *nstate, double *x_ded, double *logL, double *proba,
	     mc_law *mh, const parabox *pb, gsl_rng* rng,
	     config_base *config, FILE *FLOG, error **err);
double posterior_log_pdf_common_MPI(config_base *config, const double *x, error **err);
void run_mcmc(config_mcmc config, const parabox *pb, FILE *FLOG, FILE *CHFROM, int seed, error **err);
void mean_cov_time(double *states, int nstates, int every, config_base config, error **err);
double *burn_in_decorrelate(double *accstates, double *accparam_ded, long naccepted, long *nfinal,
			    config_mcmc config, FILE *FLOG, double **finalparam_ded, error **err);

double mean_from_chain(const double *states, long nstates, int npar, int a);
void sigma_from_chain(const double *states, int nstates, int npar, int a, double mean, double *sigma, error **err);
void mean_sigma_from_chain(config_mcmc config, const double *states, const double *states_ded, long nstates, int npar, int n_ded,
			   error **err);

double *init_param_ded(int n_ded, int nmaxchain, error **err);

void usage(int ex, const char* str, ...);

#endif

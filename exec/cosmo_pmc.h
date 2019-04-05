/* ============================================================ *
 * cosmo_pmc.h							*
 * Martin Kilbinger						*
 * ============================================================ */

#ifndef __COSMO_PMC_H
#define __COSMO_PMC_H


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sys/times.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <time.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "pmctools/errorlist.h"
#include "pmctools/mvdens.h"
#include "pmclib/tools.h"
#include "pmclib/pmc.h"
#include "pmclib/pmc_mpi.h"
#include "pmclib/parabox.h"

#include "config.h"
#include "nhist.h"
#include "param.h"
#include "stdnames.h"
#include "exec_helper.h"


#define mkpmc_base       -7000
#define mkpmc_strlen     -2 + mkpmc_base
#define mkpmc_dir        -3 + mkpmc_base

#define NSTR 128

void clean_previous_run(int iter, int niter);

void write_proposal(mix_mvdens *proposal, FILE *F1, FILE *F2, int n, const char *iterdirname, error **err);
void write_perplexity_and_ess(pmc_simu *psim, int iter, int sum_nsamples, double *sum_ess, FILE *PERP,
			      error **err);
void write_evidence(pmc_simu *psim, int iter, FILE *EVI, error **err);
void write_enc(mix_mvdens *proposal, int iter, FILE *ENC, error **err);

void histograms_and_covariance(pmc_simu *psim, const char *iterdirname, config_base *config, FILE *FLOG,
			       error **err);

void evidence_approx(config_base *config, const char *covname, const char *outname, error **err);
void revive_comp(const config_pmc *config, mix_mvdens *proposal, int i, int Nrevive, const gsl_rng *rng, error **err);
void update_and_deal_with_dead(const config_pmc *config, mix_mvdens *proposal, pmc_simu *psim,
			       int *Nrevive, const gsl_rng *rng, error **err);
void make_and_test_dir(const char *iterdirname, error **err);

void run_pmc_iteration_MPI(pmc_simu *psim, mix_mvdens **proposal_p, int iter, int this_nsamples, int *Nrevive,
			   config_pmc *config, const char *iterdirname, int myid, int nproc, int quiet,
			   gsl_rng *rng, parabox *pb, FILE *FLOG, error **err);
pmc_simu *read_pmc_iteration_MPI(mix_mvdens **proposal, int iter, int this_nsamples, config_pmc *config,
			    int myid, int nproc, FILE *PMCSIM, FILE *PROP, FILE *FLOG, error **err);
void post_processing(pmc_simu *psim, mix_mvdens *proposal, int iter, int sum_nsamples, double *sum_ess, 
		     config_base *config, const char *iterdirname, FILE *PERP, FILE *EVI, FILE *ENC,
		     FILE *FLOG, error **err);

void usage(int ex, const char* str, ...);


#endif

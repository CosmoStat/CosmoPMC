/* ============================================================ *
 * stdnames.h                                                   *
 * Martin Kilbinger 2009                                        *
 * ============================================================ */


#ifndef __STD_NAMES_H
#define __STD_NAMES_H


/* ============================================================ *
 * Common
 * ============================================================ */

/* The cosmo_pmc version.					*
 * This should be the only place in all files where the version *
 * is hard-coded.						*/
#define COSMO_PMC_VERSION "1.3"


/* Output more details to help debugging */
//#define DEBUG
#undef DEBUG

#ifdef DEBUG
#define OUTDEBUG(str) fprintf(stderr, str)
#define OUTDEBUGVA(str, ...) fprintf(stderr, str, __VA_ARGS__)
#else
#define OUTDEBUG(str)
#define OUTDEBUGVA(str, ...)
#endif

/* Log file */
#define log_name "log"

/* Random number generator initialisation */
#define rnginit_name "rng_init"

/* Likelihood histograms */
#define hist_name "chi2"


/* ============================================================ *
 * MCMC
 * ============================================================ */

/* Markov chain base name */
#define chainout "chain"

/* Chain suffixes */
#define chainall_suf "all"
#define chainacc_suf "acc"
#define chainfin_suf "fin"
#define chainpre_suf "pre"

/* (Final) chain covariance */
#define covar_name "covar"
#define inv_suf    "inv"
#define ded_suf    "+ded"

/* Directory for (updated) chain covariances */
#define covar_dir "covariances"

/* Fisher matrix */
#define fisher_name "fisher"

/* Raw Fisher matrix, no Cholesky decomposition. This matrix might not be positive or invertible! */
#define fisher_name_raw "fisher_raw"

/* Parameter mean and confidence intervals */
#define mean_name "mean"

/* Mean model parameter vector file */
#define model_mean_name "model_mean"

/* ============================================================ *
 * PMC
 * ============================================================ */

/* Directories for PMC iterations */
#define iter_dir "iter"

/* PMC sample */
#define sample_name "sample"

/* PMC simulation */
#define pmcsim_name "pmcsim"
#define pmcclip_suf "_clip"

/* Importance sampling proposal */
#define proposal_name "proposal"

/* Bayesian evience */
#define evidence_name "evidence"
#define fisher_suf      "fisher"
#define analytic_suf    "analytic"

/* Perplexity */
#define perplexity_name "perplexity"

/* Effective number of proposal components */
#define enc_name "enc"

/* Temperature */
#define temperature_name "temperature"


#endif

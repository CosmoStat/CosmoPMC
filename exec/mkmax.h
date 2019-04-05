#ifndef __MKMAX_H
#define __MKMAX_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <sys/times.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "pmctools/io.h"
#include "pmctools/maths.h"
#include "pmctools/mvdens.h"
#include "pmctools/errorlist.h"
#include "pmclib/parabox.h"
#include "pmclib/mcmc.h"
#include "pmclib/tools.h"

#include "config.h"
#include "param.h"
#include "stdnames.h"
#include "wrappers.h"
#include "all_wrappers.h"


/* Numerical derivatives: h = fh*(max-min) */
//#define fh   0.01
#define fh   0.001


double post_cg(config_base *config, const double *x, const parabox *pb, int quiet, error **err);

void dposterior_log_pdf(config_base *config, const double *x, double *df, error **err);
void mnbrak(config_base *config, double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, 
	    const double *p, const double *xi, const parabox *pb, int quiet, error **err);
double brent(config_base *config, double ax, double bx, double cx, double tol, double *xmin,
	     const double *pp, const double *xi, const parabox *pb, int quiet, error **err);
void linmin(config_base *config, double p[], double xi[], double *fret, const parabox *pb, int quiet, error **err);
void frprmn(config_base *config, double p[], double ftol, int *iter, double *fret,
	    const parabox *pb, int quiet, error **err);
double max_logP(config_base *config, double *p, const parabox *pb, double ftol, int quiet, error **err);
void set_start_parameter(double *pstate, start_t start, const double *fid, config_base config, gsl_rng *rng,
			 FILE *FLOG, error **err);
int test_maximum(const double *pstate, const config_base config, double ffh, double maxP, error **err);

double *amoeba_starting_values(double **p, int npar, config_base *config, const parabox *pb,
			       int quiet, error **err);
void amoeba(double **p, double y[], int ndim, double ftol, int *nfunk,
	    config_base *config, const parabox *pb, int quiet, error **err);
double amotry(double **p, double y[], double psum[], int ndim,
	      int ihi, double fac, config_base *config, const parabox *pb, int quiet, error **err);


#endif

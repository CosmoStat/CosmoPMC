/*
 *  nhist.h
 *  likely
 *
 *  Created by Karim Benabed on 24/02/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef __NDHIST__
#define __NDHIST__

#include "pmctools/errorlist.h"
#include "pmctools/io.h"
#include "pmctools/maths_base.h"

typedef struct {
	size_t ndim,tdim;
	double total,volume,nsamples;
  int isLog;
	size_t *nbins;
	double *limits;
	double *stps;
	double *data;
  double *var;
  void *buf;
  double *lvol;
} nd_histogram;

#define nh_base      -1800
#define nh_allocate  -1 + nh_base
#define nh_infnan    -2 + nh_base

nd_histogram  * init_nd_histogram(size_t ndim, size_t *nbins, double *limits,error **err);
nd_histogram  * init_nd_loghistogram(size_t ndim, size_t *nbins, double *limits,error **err);
void free_nd_histogram(nd_histogram **self);

void acc_histogram(size_t ndim, size_t nsamples, double * params,double *weights, size_t *pidx, nd_histogram *histo, error **err);

nd_histogram * make_histogram(size_t ndim, size_t *pidx, size_t pdim, size_t nsamples, double * params,double *weights, error **err);

double* compute_optimal_bin(size_t ndim,size_t nsamples, double *sigma, error **err);

void fill_histogram(size_t ndim, size_t nsamples, double * params,double *weights, size_t *pidx, nd_histogram *histo, error **err);

nd_histogram * build_histogram(size_t ndim, size_t *pidx, size_t pdim, size_t nsamples, double * params,double *weights,int isLog, error **err);

nd_histogram * make_kde(size_t ndim, size_t *pidx, size_t pdim, size_t nsamples, double * params,double *weights, error **err);
#endif

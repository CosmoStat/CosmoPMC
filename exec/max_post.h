#ifndef __MAX_POST_H
#define __MAX_POST_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
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

#include "nicaea/cosmo.h"
#include "nicaea/sn1a.h"

#include "param.h"
#include "pmclib/mcmc.h"
#include "pmclib/tools.h"
#include "stdnames.h"

#include "wrappers.h"
#include "all_wrappers.h"

#include "mkmax.h"
#include "exec_helper.h"


typedef enum {max_none, max_cg, max_amoeba} max_method_t;
#define smax_method_t(i) ( \
  i==max_none   ? "n" :	   \
  i==max_cg     ? "c" :	   \
  i==max_amoeba ? "a" :	   \
  "")
#define Nmax_method_t  3

void get_max_method(max_method_t *max_method, const char *optarg, error **err);
void usage(int ex, const char* str, ...);


#endif

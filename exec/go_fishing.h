#ifndef __GO_FISHING_H
#define __GO_FISHING_H

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

#include <mpi.h>

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


/* Communication tags */
#define gf_tag_size        16
#define gf_tag_fisher_part 17


/* Set DEBUG_MPI to 1 for debug outputs concerning communication between master and clients */
#define DEBUG_MPI 0
#define fprintfDEBUG(file,str,...) if (DEBUG_MPI) fprintf(file, str, __VA_ARGS__)


double fisher_element(config_base *config, const double *pos, int ab[], int quiet, error **err);
double fisher_element_adaptive(config_base *config, const double *pos, int ab[], int quiet, error **err);
mvdens *fisher_matrix(config_base *config, const double *pos, int quiet, error **err);
void fisher_matrix_part(config_base *config, const double *pos, double *fish_part,
			int my_start, int my_end, int adaptive, int quiet, error **err);
mvdens *assemble_fisher(int npar, int nproc, initial_t initial, int average_nel, int extra_nel,
			const double *fish_part_master, int my_nel_master, error **err);
void symmetrize_fisher(mvdens *fish);
void force_positive(mvdens *fish);

void usage(int ex, const char* str, ...);


#endif

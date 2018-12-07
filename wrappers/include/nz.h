/* ============================================================ *
 * nz.h								*
 * ============================================================ */


#ifndef __NZ_H
#define __NZ_H

#include "timexec.h"
#include "pmctools/maths.h"
#include "pmctools/io.h"
#include "pmctools/errorlist.h"
#include "nicaea/cosmo.h"
#include "nicaea/nofz.h"
#include "param.h"
#include "lens.h"


/* nz_fit_hist: Old CFHTLS-T0003 mode, fit explicit function n(z) to histogram data with error bars *
 * nz_read_from_files: Read n(z) data (histograms or functions) from files *
 * nz_read_from_config: Read n(z) data (histograms or functions) from config file */

typedef struct {
  int Nzbin;
  char **nzfile;
  char datname[1024], snzmode[128];
  nzmode_t nzmode;
  datcov *dc;
  char slensdata[128];
  lensdata_t lensdata;
  redshift_t *redshift;    /* To store info if read from config file */
} Nz_state;


functions_wrapper_t *init_functions_Nz(error **err);
void read_from_config_Nz(void **state, FILE *F, error **err);
void init_Nz(common_like *like, error **err);
double likeli_Nz(common_like *like, const double *param, error **err);
void print_Nz(FILE *where, void *state, error **err);

/* errors */
#define nz_base     -900
#define nz_allocate -1 + nz_base
#define nz_dead     -2 + nz_base
#define nz_choke    -3 + nz_base
#define nz_timeout  -4 + nz_base
#define nz_readerr  -5 + nz_base
#define nz_undef    -6 + nz_base


#endif /* __NZ_H */

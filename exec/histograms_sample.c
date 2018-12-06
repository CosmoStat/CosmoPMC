/* ============================================================ *
 * histograms_sample.c						*
 * Martin Kilbinger 2008-2010					*
 * Generates histograms (1d+2d) from a PMC sample.		*
 * ============================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "pmclib/pmc.h"
#include "param.h"
#include "errorlist.h"
#include "exec_helper.h"


void usage(int ex, const char* str, ...)
{
   va_list ap;

   va_start(ap, str);
   if (str!=NULL) {
      vfprintf(stderr, str, ap);
      fprintf(stderr, "\n");
   }
   va_end(ap);

   fprintf(stderr, "Usage: histograms_sample [OPTIONS] sample\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: config_pmc)\n");
   fprintf(stderr, "  -1               Only 1d histograms\n");
   fprintf(stderr, "  -2               Only 2d histograms\n");
   fprintf(stderr, "  sample           PMC sample file\n");
   fprintf(stderr, "  -h               This message\n");

   if (ex>=0) exit(ex);
}


#define NSTR 128
int main(int argc, char *argv[])
{
   FILE *F;
   char *cname, command[1024];
   error *myerr = NULL, **err;
   pmc_simu *psim;
   config_pmc config;
   double *Xall;
   int c, only_1d, only_2d;
   extern char *optarg;
   extern int optind, optopt;


   err = &myerr;


   /* Command line options */
   cname = NULL;
   only_1d = only_2d = 0;
   while ((c = getopt(argc, argv, ":c:12h")) != -1) {

      switch (c) {
	 case 'c' :
	    cname = optarg;
	    break;
	 case '1' :
	    only_1d = 1;
	    break;
	 case '2' :
	    only_2d = 1;
	    break;
	 case 'h' :
	    usage(0, NULL);
	 case ':' :
	    usage(1, "Argument for -%c option missing\n", optopt);
	 case '?' :
	    usage(2, "Unknown option -%c\n\n", optopt);
      }

   }
   if (cname==NULL) {
      cname = malloc_err(NSTR*sizeof(char), err);
      quitOnError(*err, __LINE__, stderr);
      strcpy(cname, "config_pmc");
   }
   if (argc-optind>1) usage(6, "Too many arguments");
   if (argc-optind<1) usage(7, "Sample file not given");


   /* Read config file (without initialising proposal) */
   read_config_pmc_file(&config, cname, NULL, NULL, 1, err);
   quitOnError(*err, __LINE__, stderr);

   /* Read sample and initialise PMC simulation */
   F = fopen_err(argv[optind], "r", err);
   quitOnError(*err, __LINE__, stderr);
   psim = pmc_simu_from_file(F, config.nsamples*config.fsfinal, config.base.npar, config.base.n_ded, NULL,
			     config.nclipw, err);
   quitOnError(*err, __LINE__, stderr);

   /* Merge parameters and deduced parameters */
   Xall = merge_param_param_ded(psim->X, psim->X_ded, psim->nsamples, config.base.npar,
				config.base.n_ded, err);
   quitOnError(*err, __LINE__, stderr);

   /* Mean and confidence intervals */
   mean_sigma_from_psim(psim, "mean", config.base.par, avg_mean, err);
   quitOnError(*err, __LINE__, stderr);

   /* Create histograms */
   sprintf(command, "rm -f ./chi2_*");
   system(command);
   if (! only_2d) {
      create_1d_hist(Xall, psim->nsamples, config.base.nbinhist, config.base.npar+config.base.n_ded,
		     config.base.min, config.base.max, psim->weights, ".", err);
      quitOnError(*err, __LINE__, stderr);
   }
   if (! only_1d) {
      create_2d_hist(Xall, psim->nsamples, config.base.nbinhist, config.base.npar+config.base.n_ded,
		     config.base.min, config.base.max, psim->weights, ".", err);
      quitOnError(*err, __LINE__, stderr);
   }

   covariance_from_sample_all(psim->X, psim->weights, psim->ndim, psim->nsamples,
			      Xall, psim->ndim+config.base.n_ded, ".", err);
   if (getErrorValue(*err)==mv_cholesky) {
      fprintf(stderr, "Covariance Cholesky decomposition failed, not positive, continuing...\n");
      purgeError(err);
   }
   quitOnError(*err, __LINE__, stderr);

   if (config.base.n_ded>0) free(Xall);
   pmc_simu_free(&psim);

   return 0;
}

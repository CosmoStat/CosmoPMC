/* ============================================================ *
 * meanvar_sample.c						*
 * Martin Kilbinger 2008					*
 * Calculates the mean, confidence intervals and evidence for a *
 * MCMC chain or a PMC sample (resampled simulation)		*
 * ============================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "pmclib/pmc.h"
#include "pmctools/errorlist.h"
#include "pmclib/tools.h"

#include "param.h"
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

   fprintf(stderr, "Usage: meanvar_sample [OPTIONS] sample\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: config_pmc)\n");
   fprintf(stderr, "  -w               Ignore weights (default: weights=first column of sample file)\n");
   fprintf(stderr, "  -C               Write covariance and inverse covariance to files\n");
   fprintf(stderr, "  -E               Output evidence\n");
   fprintf(stderr, "  -a AVG           Average type, AVG='mean' (default), 'median'\n");
   fprintf(stderr, "  -h               This message\n");
   fprintf(stderr, "  sample           PMC sample file\n");

   if (ex>=0) exit(ex);
}

void get_avgtype(char *savgtype, avg_t *avgtype, error **err)
{
   int i;

   STRING2ENUM(*avgtype, savgtype, avg_t, savg_t, i, Navg_t, err);
   forwardError(*err, __LINE__,);
}

#define NSTR 128
int main(int argc, char *argv[])
{
   int i, c;
   pmc_simu *psim;
   error *myerr = NULL, **err;
   FILE *PMCSIM;
   int use_weights, do_covar, do_evidence;
   char *cname, *savgtype;
   avg_t avgtype;
   config_pmc config;
   double *Xall, evi, ln_evi;
   extern char *optarg;
   extern int optind, optopt;


   err = &myerr;


   /* Command line options */
   cname       = NULL;
   savgtype    = NULL;
   use_weights = 1;
   do_covar    = do_evidence = 0;
   while ((c = getopt(argc, argv, "c:wCEa:h")) != -1) {

      switch (c) {
         case 'c' :
            cname = optarg;
            break;
         case 'w' :
            use_weights = 0;
            break;
         case 'C' :
            do_covar = 1;
            break;
         case 'E' :
            do_evidence = 1;
            break;
         case 'a' :
            savgtype = optarg;
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
   if (savgtype==NULL) {
      savgtype = malloc_err(NSTR*sizeof(char), err);
      quitOnError(*err, __LINE__, stderr);
      strcpy(savgtype, "mean");
   }
   if (argc-optind>1) usage(6, "Too many arguments");
   if (argc-optind<1) usage(7, "Sample file not given");

   /* Read config file (without initialising proposal) */
   read_config_pmc_file(&config, cname, NULL, NULL, 1, err);
   quitOnError(*err, __LINE__, stderr);

   PMCSIM      = fopen_err(argv[optind], "r", err);   quitOnError(*err, __LINE__, stderr);
   psim = pmc_simu_from_file(PMCSIM, config.nsamples*config.fsfinal, config.base.npar, config.base.n_ded, NULL,
			     config.nclipw, err);
   quitOnError(*err, __LINE__, stderr);

   if (use_weights==0) {
      for (i=0; i<psim->nsamples; i++) psim->weights[i] = 1.0/(double)psim->nsamples;
   }

   get_avgtype(savgtype, &avgtype, err);
   
   mean_sigma_from_psim(psim, "", config.base.par, avgtype, err);
   quitOnError(*err, __LINE__, stderr);

   if (do_evidence) {
      evi = evidence(psim, &ln_evi,err);
      quitOnError(*err,__LINE__,stderr);
      fprintf(stderr, "# logE lnE E\n");
      fprintf(stderr, "  % g % g % g\n", ln_evi*M_LOG10E, ln_evi, evi);
   }

   if (do_covar) {
      /* Merge parameters and deduced parameters */
      Xall = merge_param_param_ded(psim->X, psim->X_ded, psim->nsamples, config.base.npar,
				   config.base.n_ded, err);
      quitOnError(*err, __LINE__, stderr);
      covariance_from_sample_all(psim->X, psim->weights, config.base.npar, psim->nsamples, 
				 Xall, psim->ndim+config.base.n_ded, ".", err);
      quitOnError(*err, __LINE__, stderr);
   }

   pmc_simu_free(&psim);

   return 0;
}


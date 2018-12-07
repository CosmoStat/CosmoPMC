/* ============================================================ *
 * add_pmc_proposal.c						*
 * Martin Kilbinger 2008-2010					*
 *								*
 * Adds the proposal values to a PMC simulation			*
 * (files pmcsim_*). Output file is pmcsim_*+prop	        * 
 * ============================================================ */

#define Usage "add_pmc_proposal pmcsim proposal [config]\n"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "pmclib/pmc.h"
#include "pmctools/errorlist.h"
#include "pmctools/io.h"
#include "pmctools/maths.h"

#include "exec_helper.h"
#include "param.h"

void usage(int ex, const char* str, ...)
{
   va_list ap;

   va_start(ap, str);
   if (str!=NULL) {
      vfprintf(stderr, str, ap);
      fprintf(stderr, "\n");
   }
   va_end(ap);

   fprintf(stderr, "Usage: add_pmc_proposal OPTIONS FILE\n");
   fprintf(stderr, "Adds the PMC proposal value to sample points of a PMC simulation file\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: config_pmc)\n");
   fprintf(stderr, "  -p PROPOSAL      Proposal file name (required)\n");
   fprintf(stderr, "  -P               Also add posterior density\n");

   if (ex>=0) exit(ex);
}

#define NSTR 128
int main(int argc, char *argv[])
{
   FILE *F;
   config_pmc config;
   mix_mvdens *mixmvd;
   error *myerr = NULL, **err;
   char *cname, *pname;
   pmc_simu *psim;
   int i, j, c, add_post, n_ded_new;
   double logprop, *X_ded, logpost;
   extern char *optarg;
   extern int optind, optopt;


   err = &myerr;


   /* Command line options */
   cname = pname = NULL;
   add_post = 0;
   while ((c = getopt(argc, argv, ":c:p:Ph")) != -1) {

      switch (c) {
	 case 'c' :
	    cname = optarg;
	    break;
	 case 'p' :
	    pname = optarg;
	    break;
	 case 'P' :
	    add_post = 1;
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
   if (pname==NULL) usage(3, "proposal name (option -p) not given");
   if (argc-optind>1) usage(6, "Too many arguments");
   if (argc-optind<1) usage(7, "Sample file not given");


   /* Proposal (Gaussian or Student-t mixture) */
   F = fopen_err(pname, "r", err);     quitOnError(*err, __LINE__, stderr);
   mixmvd = mix_mvdens_dwnp(F, err);   quitOnError(*err, __LINE__, stderr);
   fclose(F);

   read_config_pmc_file(&config, cname, NULL, NULL, 1, err);
   quitOnError(*err, __LINE__, stderr);

   /* PMC simulation */
   F = fopen_err(argv[optind], "r", err);   
   quitOnError(*err, __LINE__, stderr);
   psim = pmc_simu_from_file(F, config.nsamples*config.fsfinal, config.base.npar, config.base.n_ded, NULL, 
			     config.nclipw, err);
   quitOnError(*err, __LINE__, stderr);
   fclose(F);

   n_ded_new = psim->n_ded + 1 + add_post;
   X_ded = malloc_err(sizeof(double)*psim->nsamples*n_ded_new, err);
   quitOnError(*err, __LINE__, stderr);

   /* Copy previous deduced parameters (if any) */
   for (i=0; i<psim->nsamples; i++) {
      for (j=0; j<psim->n_ded; j++) {
	 X_ded[i*n_ded_new+j] = psim->X_ded[i*psim->n_ded+j];
      }
   }

   /* Calculate proposal pdf */
   for (i=0; i<psim->nsamples; i++) {
      /* New deduced parameter is proposal pdf value */
      logprop = mix_mvdens_log_pdf(mixmvd, &(psim->X[i*psim->ndim]), err);
      quitOnError(*err, __LINE__, stderr);
      X_ded[i*n_ded_new+psim->n_ded] = logprop;
      if (add_post) {
	 logpost = logprop + log(psim->weights[i]);  /* log(Posterior) = log(proposal) + log(weight) */
	 //fprintf(stderr, "%g %g %g\n", prop, psim->weights[i], logpost); exit(0);
	 X_ded[i*n_ded_new+psim->n_ded+1] = logpost;
      }
   }

   psim->n_ded += n_ded_new;
   /* This creats a memory leak but X_ded cannot be freed easily because psim is one contiguous lump in memory */
   psim->X_ded = X_ded;

   sprintf(cname, "%s+prop", argv[optind]);
   psim->isLog = 0;
   out_pmc_simu_cosmo_pmc(cname, psim, config.base.par, 1.0, err);
   quitOnError(*err, __LINE__, stderr);

   pmc_simu_free(&psim);

   return 0;
}

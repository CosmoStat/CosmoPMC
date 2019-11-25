/* ============================================================ *
 * importance_sample.c						*
 * Martin Kilbinger 2009					*
 * Importance-sampling of an existing PMC simulation.		*
 * ============================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/times.h>
#include <unistd.h>

#include "pmclib/pmc.h"
#include "pmclib/pmc_mpi.h"
#include "pmclib/tools.h"
#include "pmctools/errorlist.h"

#include "param.h"
#include "exec_helper.h"

/* ============================================================ *
 * Importance sampling of the points in the PMC simulation psim *
 * with posterior as given by config.				*
 * ============================================================ */
size_t importance_sample(pmc_simu *psim, config_base *config, int quiet, error **err)
{
   int i, j, k;
   double logpost, *local_X, *local_weights, *x, *local_X_ded, MW;
   size_t nsamples_per_proc, ndim, n_ded, count;
   short *flg;
   common_like *like;

   local_X           = psim->X;
   local_X_ded       = psim->X_ded;
   local_weights     = psim->weights;
   flg               = psim->flg;
   nsamples_per_proc = psim->nsamples;
   ndim              = psim->ndim;
   n_ded             = psim->n_ded;

   MW = -1.0e30;

   for (i=count=0; i<nsamples_per_proc; i++) {
      x      = local_X+i*ndim;
      flg[i] = 0;

      /* Deal with sigma_8 */
      if (n_ded>0) {
         for (j=0; j<n_ded; j++) {
            if (config->par[ndim+j]==p_sigma8) {
               for (k=0; k<config->ndata; k++) {
                  if (config->data[k]==Lensing) {
                     like = (common_like*)(config->data_extra[k]);
                     ((lens_state*)(like->state))->sigma_8 = local_X_ded[i*n_ded+j];
                     //printf("sigma_8 = %g %p\n", ((lens_state*)(like->state))->sigma_8, &(((lens_state*)(like->state))->sigma_8));
                  }
               }
            }
         }
      }

      logpost = posterior_log_pdf_common(config, x, err);
      forwardErrorNoReturn(*err, __LINE__);
      purgeError(err);
      ParameterErrorVerb(*err, x, quiet, ndim);

      // MKDEBUG New: retrieve_ded -> &retrieve_ded according to compiler warning
      if (&retrieve_ded!=NULL && n_ded>0) {
         retrieve_ded((void*)config, &(local_X_ded[i*n_ded]), err);
         forwardErrorNoReturn(*err, __LINE__);
         purgeError(err);
         ParameterErrorVerb(*err, x, quiet, ndim);
      }

      /* Replace importance weights (log) with new log-posterior */
      local_weights[i] = logpost;

      if (i==0 || logpost>MW) {
         MW = logpost;
      }

      /* Everything fine, set flag to one and continue */
      flg[i]  = 1;
      count++;
   }

   psim->maxW = MW;

   return count;
}

void usage(int ex, const char* str, ...)
{
   va_list ap;

   va_start(ap, str);
   if (str!=NULL) {
      vfprintf(stderr, str, ap);
      fprintf(stderr, "\n");
   }
   va_end(ap);

   fprintf(stderr, "Usage: importance_sample [OPTIONS] INSAMPLE\n");
   fprintf(stderr, "Performs an importance run on a PMC sample. Run in\n");
   fprintf(stderr, " parallel with MPI (use mpirun)\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: config_pmc)\n");
   fprintf(stderr, "  -o OUTSAMPLE     Output sample name (default: 'insample.out')\n");
   fprintf(stderr, "  -q               Quiet mode\n");
   fprintf(stderr, "  -h               This message\n");
   fprintf(stderr, "  INSAMPLE         Input sample name\n");

   if (ex>=0) exit(ex);
}

#define NSTR 128
int main(int argc, char *argv[])
{
   pmc_simu *psim;
   error *myerr = NULL, **err;
   FILE *PMCSIM;
   char *cname, *inname, *outname;
   config_pmc config;
   double *weights_copy=NULL, norm, logSum_new, logSum_prev, MW;
   int i, master_samples=0, mysamples, myid, nproc, c, quiet;
   size_t nok;
   time_t t_start;
   extern char *optarg;
   extern int optind, optopt;


   err = &myerr;


   /* Command line options */
   cname = NULL;
   outname = NULL;
   quiet   = 0;
   while ((c = getopt(argc, argv, ":c:o:qh")) != -1) {

      switch (c) {
	 case 'c' :
	    cname = optarg;
	    break;
	 case 'o' :
	    outname = optarg;
	    break;
	 case 'q' :
	    quiet = 1;
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
   inname = argv[optind];
   if (outname==NULL) {
      outname = malloc_err(sizeof(char)*1024, err);  quitOnError(*err, __LINE__, stderr);
      sprintf(outname, "%s.out", inname);
   }


   /* Initialise  MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
   if (myid==0) {
      fprintf(stderr, "%s compiled on %s %s\n", __FILE__, __DATE__, __TIME__);
      time(&t_start);
      fprintf(stderr, "started at %s\n", ctime(&t_start));
      if (! quiet) fprintf(stderr, "Number of nodes nproc = %d\n", nproc);
   }


   /* Read config file */
   sleep(myid*0.01);
   read_config_pmc_file(&config, cname, NULL, NULL, 0, err);
   quitOnError(*err, __LINE__, stderr);


   if (myid==0) {        /* Master */

      /* Read PMC simulation */
      PMCSIM      = fopen_err(inname, "r", err);
      quitOnError(*err, __LINE__, stderr);
      psim = pmc_simu_from_file(PMCSIM, config.nsamples*config.fsfinal, config.base.npar, config.base.n_ded, NULL,
				config.nclipw, err);
      quitOnError(*err, __LINE__, stderr);
      fclose(PMCSIM);

      /* Now, the weights are normalized and log-weights */

      logSum_prev = psim->logSum;
      testErrorRetVA(psim->ndim != config.base.npar, mk_npar,
		     "Number of parameters not consistent between config file (%s, npar=%d) and "
		     "pmc sim file (%s, npar=%d)", *err, __LINE__, mk_npar, cname, config.base.npar,
		     inname, psim->ndim);

      /* Copy previous log-weights */
      weights_copy = malloc_err(sizeof(double)*psim->nsamples, err);
      quitOnError(*err, __LINE__, stderr);
      memcpy(weights_copy, psim->weights, sizeof(double)*psim->nsamples);

   } else {            /* Clients */

      /* Init psim container with dummy size */
      psim = pmc_simu_init_plus_ded(1, config.base.npar, config.base.n_ded, err);
      quitOnError(*err, __LINE__, stderr);

   }


   /* Send and receive simulations */
   if (myid==0) {    /* Master */

      if (! quiet) fprintf(stderr, "proc %d >> send pmc sim\n", myid);
      mysamples      = send_simulation(psim, nproc, err);              quitOnError(*err, __LINE__, stderr);
      master_samples = psim->nsamples;
      psim->nsamples = mysamples;

   } else {          /* Clients */

      if (! quiet) fprintf(stderr, "proc %d >> receive pmc sim\n", myid);
      receive_simulation(psim, nproc, myid, err);
      quitOnError(*err, __LINE__, stderr);

   }


   /* Loop over all simulated points (master and clients) */
   if (! quiet) fprintf(stderr, "proc %d >> work on %ld samples\n", myid, psim->nsamples);
   nok = importance_sample(psim, &(config.base), quiet, err);
   quitOnError(*err, __LINE__, stderr);
   if (! quiet) fprintf(stderr, "proc %d >> finished importance weights, nok=%zu\n", myid, nok);

   /* Send and receive importance weights */
   if (myid!=0) {    /* clients */

      /* Send importance weights */
      if (! quiet) fprintf(stderr, "proc %d >> send weights\n" ,myid);
      send_importance_weight(myid, nproc, psim, nok);

   } else {     /* master */

      /* Receive computation from clients */
      mysamples      = psim->nsamples;
      psim->nsamples = master_samples;     
      if (! quiet) fprintf(stderr, "proc %d >> receive weights\n", myid);
      nok = receive_importance_weight(psim, nproc, nok, mysamples, err);
      quitOnError(*err, __LINE__, stderr);

      /* MKDEDEBUG New: Normalize (which creates weights), take log again */
      psim->isLog = 1;
      norm = normalize_importance_weight(psim, err);
      quitOnError(*err, __LINE__, stderr);
      logSum_new = psim->logSum;
      for (i=0; i<psim->nsamples; i++) {
	 psim->weights[i] = log(psim->weights[i]);
      }

      if (! quiet) fprintf(stderr, "proc %d >> nsamples = %ld\n", myid, psim->nsamples);

      /* Add previous (copied above) log-weights to new (replaced in importance_sample) log-weights *
       * Since weights are normalized, add both log-sums to weights.                                */
      MW = psim->maxW;
      for (i=0; i<psim->nsamples; i++) {
	 psim->weights[i] += log(weights_copy[i]) + logSum_prev + logSum_new;

	 /* New in v > 1.2: Update maximum weight */
	 if (psim->weights[i] > MW) MW = psim->weights[i];
      }
      psim->maxW = MW;

      /* Normalise weights */
      if (! quiet) fprintf(stderr, "Master: normalize_importance_weight\n");
      psim->isLog = 1;
      norm = normalize_importance_weight(psim, err);
      quitOnError(*err, __LINE__, stderr);

      /* Print PMC simulation */
      out_pmc_simu_cosmo_pmc(outname, psim, config.base.par, norm, err);
      quitOnError(*err, __LINE__, stderr);

      free(weights_copy);

   }

   pmc_simu_free(&psim);


   if (myid==0) {       /* master */
      end_time(t_start, stderr);
      fprintf(stderr, "importance_sample finished\n");
   }

   MPI_Finalize();
 
   return 0;
}

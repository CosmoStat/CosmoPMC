/* ============================================================ *
 * Martin Kilbinger 2008					*
 * Reads a sample (MCM chain or PMC (re)sample and prints the   *
 * corresponding w(theta) with error for each parameter vector.	*
 * ============================================================ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/times.h>
#include <time.h>
#include <math.h>
#include <sys/times.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "errorlist.h"
#include "param.h"
#include "config.h"

#include "nicaea/hod.h"
#include "halo.h"

#include "pmclib/pmc.h"
#include "exec_helper.h"


int same_parameter(const double *X1, const double *X2, int npar)
{
   int i;

   for (i=0; i<npar; i++) {
      if (fabs(X1[i]-X2[i])>EPSILON1) return 0;
   }
   return 1;
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

   fprintf(stderr, "Usage: sample_GalCorr [OPTIONS] sample\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: config_pmc)\n");
   fprintf(stderr, "  sample           PMC sample file\n");

   if (ex>=0) exit(ex);
}

void get_w_theta(cosmo_hm *def, const double *params, const par_t *par, int npar, double *theta, size_t Ntheta, double *wth, error **err)
{
   int j, i_bin, j_bin;
   cosmo_hm *halomodel;
   double *w1hgcs, *w1hgss, *w2hg;

   halomodel  = copy_parameters_hm(def, err);                       forwardError(*err, __LINE__,);
   fill_parameters_hm(halomodel, params, par, npar, err);           forwardError(*err, __LINE__,);
   updateFrom_hm(def, halomodel, err);                              forwardError(*err, __LINE__,);

   i_bin = j_bin = 0;
   w1hgcs = woftheta(halomodel, p1hgcs, theta, Ntheta, i_bin, j_bin, err); forwardError(*err, __LINE__,);
   w1hgss = woftheta(halomodel, p1hgss, theta, Ntheta, i_bin, j_bin, err); forwardError(*err, __LINE__,);
   w2hg   = woftheta(halomodel, p2hg, theta, Ntheta, i_bin, j_bin, err);   forwardError(*err, __LINE__,);

   for (j=0; j<Ntheta; j++) {
      wth[j]  = w1hgcs[j] + w1hgss[j] + w2hg[j];
   }

   free(w1hgcs);
   free(w1hgss);
   free(w2hg);
   free_parameters_hm(&halomodel);
}

void get_NofM(cosmo_hm *def, const double *params, const par_t *par, int npar, const double *M, size_t NlogM, double *NofM, error **err)
{
   int j;
   cosmo_hm *halomodel;

   halomodel  = copy_parameters_hm(def, err);                       forwardError(*err, __LINE__,);
   fill_parameters_hm(halomodel, params, par, npar, err);           forwardError(*err, __LINE__,);
   updateFrom_hm(def, halomodel, err);                              forwardError(*err, __LINE__,);

   for (j=0; j<NlogM; j++) {
     NofM[j] = Ngal(halomodel, M[j], def->log10Mstar_min, def->log10Mstar_max, err); forwardError(*err, __LINE__,);
   }

   free_parameters_hm(&halomodel);
}


void write_mean_and_var(const char *name, const double *x, double *y, size_t N, const pmc_simu *psim, error **err)
{
   const double conf_123_half[3] = {conf_68/2.0, conf_95/2.0, conf_99/2.0};

   double *pmean, *pvar;
   int i, j;
   FILE *F;

   pmean = malloc_err(sizeof(double)*N, err);          forwardError(*err, __LINE__,);
   pvar  = calloc_err(N*6, sizeof(double), err);       forwardError(*err, __LINE__,);

   F = fopen_err(name, "w", err);                      forwardError(*err, __LINE__,);
   for (j=0; j<N; j++) {

      pmean[j] = mean_from_psim(y, psim->weights, psim->flg, psim->nsamples, N, j);
      sigma_from_psim(y, psim->weights, psim->flg, psim->nsamples, N, j, pmean[j], pvar+6*j, conf_123_half, err);
      forwardError(*err, __LINE__,);

      fprintf(F, "%g  % g", x[j], pmean[j]);
      for (i=0; i<3; i++) fprintf(F, "   % 9.5g % 9.5g", pvar[6*j+i], pvar[6*j+i+3]);
      fprintf(F, "\n");

   }
   fclose(F);

   free(pmean);
   free(pvar);
}

#define NSTR 128
int main(int argc, char *argv[])
{

   config_pmc config;
   error *myerr = NULL, **err;
   parabox *pb;
   pmc_simu *psim;
   int i, j;
   FILE *F;
   common_like *like;
   halo_state *halostate;
   double *wth, *NofM, logth_min, logth_max, delta, *theta, *M, logM_min, logM_max, delta_logM;
   char *cname;
   size_t Ntheta, NlogM;
   int c;
   extern char *optarg;
   extern int optind, optopt;


   err = &myerr;


   /* Command line options */
   cname = NULL;
   while ((c = getopt(argc, argv, ":c:h")) != -1) {

      switch (c) {
	 case 'c' :
	    cname = optarg;
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


   /* Reading header of config file */
   read_config_pmc_file(&config, cname, NULL, NULL, 1, err);
   quitOnError(*err, __LINE__, stderr);

   /* Initialise parameter box */
   pb = parabox_from_config(config.base.npar, config.base.min, config.base.max, err);
   quitOnError(*err, __LINE__, stderr);

   /* Initialise PMC simulation */
   F        = fopen_err(argv[optind], "r", err);
   quitOnError(*err, __LINE__, stderr);
   /* NEW: with clip weights */
   psim = pmc_simu_from_file(F, config.nsamples, config.base.npar, config.base.n_ded, NULL, config.nclipw, err);
   quitOnError(*err, __LINE__, stderr);
   fclose(F);

   /* Default halomodel from file */
   like = (config.base.data_extra[0]);
   halostate = (halo_state*)like->state;

   if (strcmp(halostate->model_file, "-")!=0) {
      F = fopen_err(halostate->model_file, "r", err);
      quitOnError(*err, __LINE__, stderr);
      read_cosmological_parameters_hm(&halostate->model, F, err);
      quitOnError(*err, __LINE__, stderr);
      fclose(F);
   } else {
      halostate->model = set_cosmological_parameters_to_default_hm(err);
      quitOnError(*err, __LINE__, stderr);
   }

   Ntheta    = 60;
   logth_min = log10(0.0001);
   logth_max = log10(10.0);

   wth       = malloc_err(sizeof(double) * psim->nsamples * Ntheta, err); quitOnError(*err, __LINE__, stderr);
   delta     = (logth_max - logth_min) / (double)Ntheta;
   theta     = malloc_err(sizeof(double) * Ntheta, err);                  quitOnError(*err, __LINE__, stderr);
   for (j=0; j<Ntheta; j++) {
      theta[j] = pow(10.0, logth_min + (double)j * delta);
   }

   NlogM      = 30;
   logM_min   = 10.0;   /* log10 */
   logM_max   = 16.0;

   NofM       = malloc_err(sizeof(double) * psim->nsamples * NlogM, err);  quitOnError(*err, __LINE__, stderr);
   delta_logM = (logM_max - logM_min) / (double)NlogM;
   M          = malloc_err(sizeof(double) * NlogM, err);  quitOnError(*err, __LINE__, stderr);
   for (j=0; j<NlogM; j++) {
      M[j] = pow(10.0, logM_min + (double)j * delta_logM);
   }

   for (i=0; i<psim->nsamples; i++) {

      printf("%d ", i); fflush(stdout);
      get_w_theta(halostate->model, psim->X + i*config.base.npar, config.base.par, config.base.npar, theta, Ntheta, wth + i*Ntheta, err);
      quitOnError(*err, __LINE__, stderr);

      get_NofM(halostate->model, psim->X + i*config.base.npar, config.base.par, config.base.npar, M, NlogM, NofM + i*NlogM, err);
      quitOnError(*err, __LINE__, stderr);

   }
   printf("\n");

   write_mean_and_var("w_theta.mean", theta, wth, Ntheta, psim, err); quitOnError(*err, __LINE__, stderr);
   write_mean_and_var("NofM.mean", M, NofM, NlogM, psim, err); quitOnError(*err, __LINE__, stderr);


   free(theta);
   free(wth);
   free(M);
   free(NofM);

   pmc_simu_free(&psim);
   free_parabox(&pb);

   return 0;
}
#undef NSTR


/* ============================================================ *
 * exec_helper.c						*
 * Martin Kilbinger, Karim Benabed				*
 * ============================================================ */


#include "exec_helper.h"


/* ============================================================ *
 * Returns GSL random number generator. Initialised with seed.  *
 * If seed==-1, initialised with current time.			*
 * ============================================================ */
gsl_rng *init_random(int seed, FILE *FLOG)
{
   gsl_rng *rng;
   struct tms Time;
   clock_t cl;
   unsigned u;

   if (seed==-1) {
      cl = times(&Time);
      u  = (unsigned)cl;
   } else {
      u  = (unsigned)seed;
   }

   rng = gsl_rng_alloc(gsl_rng_default);
   gsl_rng_set(rng, u);

   if (FLOG) {
      fprintf(FLOG, "Seed for random generator = %u\n", u);
   }

   return rng;
}

void print_mean_sigma(const char *pname, const double *pmean, const double *pvar, const par_t *par, int ntot,
                      avg_t avgtype, error **err)
{
   int len, i, j;
   FILE *F;

   len = strlen(pname);
   if (len==0) {
      F = stdout;
   } else {
      F = fopen_err(pname, "w", err);
      forwardError(*err,__LINE__,);
   }
   fprintf(F, "# i   spar                % 6s +68.27/2%% -68.27/2%%   +95.45/2%% -95.45/2%%   +99.73/2%% -99.73/2%%\n",
            savg_t(avgtype));
   for (i=0; i<ntot; i++) {
      fprintf(F, " %2d   %-15s %10.5g", i, spar_t(par[i]), pmean[i]);
      for (j=0; j<3; j++) fprintf(F, "   % 9.5g % 9.5g", pvar[6*i+j], pvar[6*i+j+3]);
      fprintf(F, "\n");
   }
   if (len!=0) {
      fclose(F);
   }
}

double *mean_sigma_from_psim(pmc_simu *psim, const char *pname, const par_t *par, avg_t avgtype, error **err)
{
  const double conf_123_half[3] = {conf_68/2.0, conf_95/2.0, conf_99/2.0};
  double *pmean, *pvar;
  int i, ntot;
  
  testErrorRet(psim->isLog==1, pmc_isLog, "Weights are in log, transform them before calculating mean and errors",
	       *err, __LINE__, NULL);
  
  ntot  = psim->ndim+psim->n_ded;
  pmean = calloc_err(ntot, sizeof(double), err);     forwardError(*err, __LINE__, NULL);
  pvar  = calloc_err(ntot*6, sizeof(double), err);   forwardError(*err, __LINE__, NULL);

  for (i=0; i<psim->ndim; i++) {

     switch (avgtype) {
        case avg_mean :
           pmean[i] = mean_from_psim(psim->X, psim->weights, psim->flg, psim->nsamples, psim->ndim, i);
           break;
        case avg_median :
           pmean[i] = median_from_psim(psim->X, psim->weights, psim->flg, psim->nsamples, psim->ndim, i, err);
           forwardError(*err, __LINE__, NULL);
           break;
     }

     sigma_from_psim(psim->X, psim->weights, psim->flg, psim->nsamples,
     		     psim->ndim, i, pmean[i], pvar+6*i, conf_123_half, err);
     forwardError(*err, __LINE__, NULL);
  }

  /* Deduced parameters */
  for (i=0; i<psim->n_ded; i++) {

     switch (avgtype) { 
        case avg_mean :
           pmean[psim->ndim+i] = mean_from_psim(psim->X_ded, psim->weights, psim->flg, psim->nsamples,
                 psim->n_ded, i);
           break;
        case avg_median :
           pmean[psim->ndim+i] = median_from_psim(psim->X_ded, psim->weights, psim->flg, psim->nsamples,
                 psim->n_ded, i, err);
           forwardError(*err, __LINE__, NULL);
           break;
     }

     sigma_from_psim(psim->X_ded, psim->weights, psim->flg, psim->nsamples,
           psim->n_ded, i, pmean[psim->ndim+i], pvar+6*(i+psim->ndim), conf_123_half, err);
     forwardError(*err, __LINE__, NULL);
  }

  print_mean_sigma(pname, pmean, pvar, par, ntot, avgtype, err);
  forwardError(*err, __LINE__, NULL);

  free(pvar);

  return pmean;
}


void write_mean(double *pmean, const char *basename, int ndata, void **data_extra, error **err)
{
   int i;
   double logL;
   common_like *like;
   char name[128];

   for (i=0; i<ndata; i++) {
      like = (common_like*)(data_extra[i]);
      like->want_model = 1;
      logL = like->func_wrapper->func_likeli(like, pmean, err);
      forwardError(*err, __LINE__,);

      sprintf(name, "%s_%d", basename, i);
      write_model_raw(like->model, like->Nmodel, name, err);
   }
}

void write_model_raw(double *model, int Nmodel, const char *name, error **err)
{
   FILE *F;
   int i;

   F = fopen_err(name, "w", err);    forwardError(*err, __LINE__,);
   for (i=0; i<Nmodel; i++) {
      fprintf(F, "%d %g\n", i, model[i]);
   }
   fclose(F);
}


int par1dw_cmp(const void *av, const void *bv)
{
  const struct par1dw *a, *b;
  a = (const struct par1dw*)av;
  b = (const struct par1dw*)bv;
  if (a->par<b->par) return -1;
  if (a->par>b->par) return +1;
  return 0;
}

/* Returns median of a-th parameter */
double median_from_psim(double *X, double *w, short *flg, int nsamples, int ndim, int a, error **err)
{
  double m;
  int i, n;
  struct par1dw *tmp;


  /* Copy relevant data and sort according to parameter values */
  tmp = malloc_err(sizeof(struct par1dw)*nsamples, err);
  forwardError(*err, __LINE__, 0.0);

  for (i=0,n=0; i<nsamples; i++) {
    if (flg[i]==0) continue;
    tmp[i].par = X[i*ndim+a];
    tmp[i].w   = w[i];
    n++;
  }
  qsort(tmp, n, sizeof(struct par1dw), par1dw_cmp);

  /* TODO(?): check isLog */
  for (i=0,m=0.0; i<nsamples; i++) {
     m += tmp[i].w;

     /* This case should not occur due to numerical inaccuracies */
     if (m==0.5) return tmp[i].par;

     /* Beyond half density: return mean of current and previous point */
     if (m>0.5) {
        return (tmp[i].par + tmp[i-1].par)/2;
     }
  }

  *err = addErrorVA(err_median, "Sum of normalised weights does not go beyond 0.5",
         *err, __LINE__, 0.0);
   return 0;
}

void sigma_from_psim(double *X, double *weights, short *flg, int nsamples, int ndim,
      int a, double mean, double *sigma, const double confidence[], error **err)
{
  /* erf({1,2,3}/sqrt(2)) */
  int i, j, imean, n;
  double sum;
  struct par1dw *tmp;
  
  /* Copy relevant data and sort according to parameter values */
  tmp = malloc_err(sizeof(struct par1dw)*nsamples, err);
  forwardError(*err, __LINE__,);
  
  for (i=0,n=0; i<nsamples; i++) {
    if (flg[i]==0) continue;
    tmp[i].par = X[i*ndim+a];
    tmp[i].w   = weights[i];
    n++;
  }
  qsort(tmp, n, sizeof(struct par1dw), par1dw_cmp);
  
  /* Look for mean */
  i = 0;
  while (tmp[i].par<mean) {
    i++;
    if (i==n-1) break;
  }
  imean = i;
  
  /* Sum up right of mean until confidence volume */
  for (j=0; j<3; j++) {
    sum = 0.0;
    i = imean;
    while (sum<=confidence[j] && i<n) {
      sum += tmp[i].w;
      i++;
    }
    sigma[j] = tmp[i].par-mean;
    /* Negative value: boundary hit before cl volume reached */
    if (i==n) sigma[j] = -1.0;
  }

  /* MKDEBUG */
  /* printf("MKDEBUG par=%d mean=%g i=%d\n", a, tmp[imean].par, imean);
  sum = 0.0;
  i = imean;
  while (i<n) {
   sum += tmp[i].w;
   i++;
  }
  printf("MKDEBUG par=%d, sum right = %g\n", a, sum); */

  /* Sum up left of mean until confidence volume */
  for (j=0; j<3; j++) {
    sum = 0.0;
    i = imean-1;
    while (sum<=confidence[j] && i>=0) {
      sum += tmp[i].w;
      //if (a==0) printf("MKDEBUG par=%d cf=%d %g<=%g  %g %g\n", a, j, sum, confidence[j], tmp[i].par, tmp[i].w);
      i--;
    }
    sigma[j+3] = mean-tmp[i].par;
    if (i==-1) sigma[j+3] = -1.0;
  }

  /* MKDEBUG */
  /* sum = 0.0;
  i = imean-1;
  while (i>=0) {
   sum += tmp[i].w;
   i--;
  }
  printf("MKDEBUG par=%d, sum left = %g\n", a, sum); */

  free(tmp);
}


void confidence_level(double *X, double *weights, short *flg, int nsamples,
                      int ndim, int a, double val, int dir, const double *conf_xxx, double *cl, error **err)
{
  int i;
  double pvar[6];
  
  testErrorRet(abs(dir)!=1, pmc_undef, "Direction has to be +1 or -1", *err, __LINE__,);
  
  /* Call sigma_from_psim with val instead of pmean */
  sigma_from_psim(X, weights, flg, nsamples, ndim, a, val, pvar, conf_xxx, err);
  
  for (i=0; i<3; i++) {
    cl[i] = pvar[i + (1-dir)/2*3];
  }
} 

/* ============================================================ *
 * Calculates and writes to files the covariance and inverse    *
 * of all sample and sample+deduced parameters.			*
 * ============================================================ */
#define NSTR 128
void covariance_from_sample_all(const double *X, const double *weights, int ndim, long nsamples,
                                const double *X_plus_ded, int ndim_plus_ded, const char *dir, error **err)
{
  char name[NSTR], name_inv[NSTR];
  
  sprintf(name, "%s/%s.%s", dir, covar_name, chainfin_suf);
  sprintf(name_inv, "%s/%s%s.%s", dir, covar_name, inv_suf, chainfin_suf);
  covariance_from_sample(X, weights, ndim, nsamples, name, name_inv, err);
  forwardError(*err, __LINE__,);
  
  sprintf(name, "%s/%s%s.%s", dir, covar_name, ded_suf, chainfin_suf);
  sprintf(name_inv, "%s/%s%s%s.%s", dir, covar_name, inv_suf, ded_suf, chainfin_suf);
  covariance_from_sample(X_plus_ded, weights, ndim_plus_ded, nsamples, name, name_inv, err);
  forwardError(*err, __LINE__,);
}
#undef NSTR

/* ============================================================ *
 * Front-end for estimate_param_covar_weight. Writes covariance *
 * and inverse to files.					*
 * ============================================================ */
void covariance_from_sample(const double *X, const double *weight, int ndim, long nsample,
                            const char *name, const char *cinvname, error **err)
{
  double *pmean, *pvar;		   /* Mean and variance of final sample */
  mvdens *cov;		 	   /* Covariance of final sample (unweighted) */
  FILE *F;

  pmean   = calloc_err(ndim, sizeof(double), err);
  forwardError(*err, __LINE__, );
  pvar    = calloc_err(ndim*ndim, sizeof(double), err);
  forwardError(*err, __LINE__,);

  estimate_param_covar_weight(ndim, nsample, 0, X, weight, pmean, pvar, err);
  forwardError(*err, __LINE__,);

  cov = mvdens_alloc(ndim, err);         forwardError(*err, __LINE__,);
  mvdens_from_meanvar(cov, pmean, pvar, 1.0);
  free(pmean);
  free(pvar);

  F = fopen_err(name, "w",err); forwardError(*err,__LINE__,);
  mvdens_dump(F, cov);
  fclose(F);

  mvdens_inverse(cov, err);               forwardError(*err, __LINE__,);
  //sm2_inverse(cov->std, ndim, err);       forwardError(*err, __LINE__,);
  F = fopen_err(cinvname, "w",err);       forwardError(*err,__LINE__,);
  mvdens_dump(F, cov);
  fclose(F);
}

void write_header_cosmo_pmc(FILE *OUT, const par_t *par, int npar, int n_ded)
{
   int i;

   fprintf(OUT, "# npar = %d, n_ded = %d\n", npar, n_ded);

   //fprintf(OUT, "# number   1 %d         param\n", npar);
   //fprintf(OUT, "# number   1 %d         param_ded\n", n_ded);
   fprintf(OUT, "#         weight            chi2");
   for (i=0; i<npar; i++) {
      fprintf(OUT, "        param[%d]", i);
   }
   for (i=0; i<n_ded; i++) {
      fprintf(OUT, "    param_ded[%d]", i);
   }
   fprintf(OUT, "\n");

   fprintf(OUT, "#                               ");
   for (i=0; i<npar; i++) {
      fprintf(OUT, "% 16s", spar_t(par[i]));
   }
   for (i=0; i<n_ded; i++) {
      //fprintf(OUT, "    param_ded[%d]", i);
   }
   fprintf(OUT, "\n");
}

void print_step_cosmo_pmc(FILE* where, double weight, double loglkl, size_t npar, size_t n_ded, const double *params,
		   const double *params_ded)
{
   FILE *rf;
   size_t i;
	    
   rf = where;
   if (where==NULL) rf = stdout;
   /* MCMC: accept -logL
      PMC:  weight -component */
   fprintf(rf, "%16.9g%16.9g", weight, -loglkl);
   for (i=0;i<npar;i++) {
      fprintf(rf, "%16.9g", params[i]);
   }
   for (i=0;i<n_ded;i++) {
      fprintf(rf, "%16.9g", params_ded[i]);
   }
   fprintf(rf, "\n"); fflush(rf);
}

/* Prints PMC simulation in CosmoPMC format */
void out_pmc_simu_cosmo_pmc(const char *name, const pmc_simu *psim, const par_t *par, double norm, error **err)
{
   FILE *OUT;
   int i;
   double logw, lognorm;

   OUT = fopen_err(name, "w",err);
   forwardError(*err,__LINE__,);

   lognorm = psim->logSum;

   write_header_cosmo_pmc(OUT, par, psim->ndim, psim->n_ded);
   for (i=0; i<psim->nsamples; i++) {

      if (psim->flg[i]==0) continue;

      logw = psim->weights[i];
      if (!psim->isLog) logw = log(logw);


      /* Back to unnormalised weights (log) */
      logw += lognorm;

      print_step_cosmo_pmc(OUT, logw, psim->indices[i], psim->ndim, psim->n_ded, &(psim->X[i*psim->ndim]), 
			   &(psim->X_ded[i*psim->n_ded]));

   }
   fclose(OUT);
}


void make_and_test_dir(const char *dirname, error **err)
{
   int nok;

   if(strcmp(dirname,".")){
     mkdir(dirname, 0755);
     nok = chdir(dirname);
     if (nok==0) {    /* chdir successful */ 
       chdir("..");
     } else {
       *err = addErrorVA(err_dirname, "Could not create and change into directory %s. Make sure you have "
			 "write permission.", *err, __LINE__, dirname);
       return;
     }
   }
   
}

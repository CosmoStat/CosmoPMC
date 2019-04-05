/* ============================================================ *
 * param.c							*
 * Martin Kilbinger 2007-2009					*
 * ============================================================ */


#include "param.h"

/* In case of error: Prints error to files F1 and F2 (can be NULL), purges it and returns 1. *
 * Else: No action, returns 0.								     *
 * params can be NULL.									     */
int checkPrintErr_and_continue(error **err, FILE *F1, FILE *F2, const double *params, int npar)
{
   if (isError(*err)) {
      if (F1!=NULL) {
         fprintf(F1, "*** Error:\n");
         if (params!=NULL) print_parameter(F1, npar, params);
         printError(F1, *err);
         fprintf(F1, "continuing...\n");
         fflush(F1);
      }	    

      if (F2!=NULL) {
         fprintf(F2, "*** Error:\n");
         if (params!=NULL) print_parameter(F2, npar, params);
         printError(F2, *err);
         fprintf(F2, "continuing...\n");
         fflush(F2);
      }

      purgeError(err);
      return 1;
   }

   return 0;
}

/* ============================================================ *
 * Reads the prior and returns it if given in the config file.	*
 * Otherwise returns NULL.					*
 * ============================================================ */
mvdens *assign_prior(char sprior[], int npar, int no_init, error **err)
{
   FILE *PRIOR;
   mvdens *prior;

   if (strcmp(sprior, "-")!=0) {

      if (no_init == 1) {
         /* Dummy prior */
         prior = mvdens_alloc(1, err);         forwardError(*err, __LINE__, NULL);
      } else {
         PRIOR = fopen_err(sprior, "r", err);     forwardError(*err, __LINE__, NULL);
         prior = mvdens_dwnp(PRIOR, err);         forwardError(*err, __LINE__, NULL);
         fclose(PRIOR);
         mvdens_cholesky_decomp(prior, err);      forwardError(*err, __LINE__, NULL);
      }

      //testErrorRetVA(npar!=prior->ndim, mk_npar, "Wrong dimension %d of prior (file '%s'), expected %d",
      //	   *err, __LINE__, NULL, prior->ndim, sprior, npar);
   } else {

      prior = NULL;

   }

   return prior;
}

/* ============================================================ *
 * Reads the base part of a config file (common to pmc, mcmc).  *
 * ============================================================ */
void read_config_base(config_base *config, FILE *F, int no_init, error **err)
{
   config_element c = {0, 0.0, ""};
   int i, j, n, n_tot;
   char s[500];

   /* ============================================================ *
    * Version							   *
    * ============================================================ */

   CONFIG_READ(config, version, d, F, c, err);

   if (config->version > atof(COSMO_PMC_VERSION)) {
      /*
      fprintf(stderr, "Warning: Config file version (%g) is higher than code version (%g).\n"
	      "Upwards compatibility of the code cannot be guaranteed, continue at own risk\n",
	      config->version, atof(COSMO_PMC_VERSION));
      */
   } else if (config->version < atof(COSMO_PMC_VERSION)) {
      /* For future versions: Check whether code is downwards compatible */
   }


   /* ============================================================ *
    * Parameter section						   *
    * ============================================================ */

   /* Number of parameters, deduced parameters */
   CONFIG_READ(config, npar, i, F, c, err);
   CONFIG_READ(config, n_ded, i, F, c, err);

   n_tot = config->npar + config->n_ded;

   config->spar = malloc_err(sizeof(char*)*n_tot, err);
   forwardError(*err, __LINE__,);
   for (i=0; i<n_tot; i++) {
      config->spar[i] = malloc_err(sizeof(char)*50, err);                 forwardError(*err, __LINE__,);
   }
   CONFIG_READ_S_ARR(config, spar, s, i, n_tot, F, c, err);

   config->par = malloc_err(sizeof(par_t)*n_tot, err);                    forwardError(*err, __LINE__,);
   for (i=0; i<n_tot; i++) {
      STRING2ENUM(config->par[i], config->spar[i], par_t, spar_t, j, Npar_t, err);
   }

   config->min = malloc_err(sizeof(double)*n_tot, err);                   forwardError(*err, __LINE__,);
   config->max = malloc_err(sizeof(double)*n_tot, err);                   forwardError(*err, __LINE__,);

   CONFIG_READ_ARR(config, min, d, i, n_tot, s, F, c, err);
   CONFIG_READ_ARR(config, max, d, i, n_tot, s, F, c, err);

   for (j=0,config->logpr_default=0.0; j<config->npar; j++) {
      testErrorRetVA(config->max[j]<=config->min[j], mcmc_infnan,
		     "The minumum #%d (%g) is not smaller than the maximum (%g), check config file",
		     *err, __LINE__,, j, config->min[j], config->max[j]);
      config->logpr_default -= log(config->max[j] - config->min[j]);
   }
   testErrorRet(!finite(config->logpr_default), mcmc_infnan,
		"logpr_default not finite. Check min, max in config file", *err, __LINE__,);


   /* ============================================================ *
    * Data section						   *
    * ============================================================ */

   /* Number and type of data sets */
   CONFIG_READ(config, ndata, i, F, c, err);

   testErrorRet(config->ndata<=0, mcmc_negative, "ndata (in config file) must be larger than zero",
		*err, __LINE__,);
   config->data = (data_t*)malloc_err(sizeof(data_t)*config->ndata, err);
   forwardError(*err, __LINE__,);
   for (i=0; i<config->ndata; i++) {
      config->data[i] = -1;
      CONFIG_READ_S(config, sdata, s, F, c, err);
      STRING2ENUM(config->data[i], config->sdata, data_t, sdata_t, j, Ndata_t, err);
   }

   /* Extra for each dataset */
   config->data_extra = (common_like**)malloc_err(sizeof(common_like*)*config->ndata, err);
   forwardError(*err, __LINE__,);
   for (i=0; i<config->ndata; i++) {
      config->data_extra[i] = extra_init_common(config->data[i], config->par, config->npar, err);
      forwardError(*err, __LINE__,);

      read_extra_for_data(F, config->data[i], config->data_extra[i], no_init, err);
      forwardError(*err, __LINE__,);
   }


   /* Prior */
   CONFIG_READ_S(config, sprior, s, F, c, err);
   config->prior = assign_prior(config->sprior, config->npar, no_init, err);      forwardError(*err, __LINE__,);
   if (config->prior!=NULL) {
      /* Number of parameters to which priors applies */
      CONFIG_READ(config, nprior, i, F, c, err);
      /*
      testErrorRetVA(no_init != 1 && config->nprior!=config->prior->ndim, mk_npar,
		     "Inconsistent dimension %d of prior (file '%s') and config-file key 'nprior' %d",
		     *err, __LINE__,, config->prior->ndim, config->sprior, config->nprior);
		     */

      /* Flag array for prior variables */
      config->indprior = malloc_err(sizeof(int)*config->npar, err);     forwardError(*err, __LINE__,);
      CONFIG_READ_ARR(config, indprior, i, i, config->npar, s, F, c, err);
      for (i=n=0; i<config->npar; i++) {
	 n += config->indprior[i];
	 testErrorRetVA(config->indprior[i]!=0 && config->indprior[i]!=1, mcmc_prior,
			"Prior flag array (indprior in config file) has to be 0 or 1,\n"
			"%d found for entry #%d", *err, __LINE__,, config->indprior[i], i);
      }

      testErrorRetVA(n!=config->nprior, mcmc_dimension,
		   "Number of entries in prior flag array (%d) not consistent with first entry (%d)\n"
		     "(sprior in config file)", *err, __LINE__,, n, config->nprior);

      if (no_init == 1) {
         free(config->prior);
         config->prior = NULL;
         free(config->indprior);
         config->indprior = NULL;
      }

   } else {
      config->indprior = NULL;
   }

}

/* ============================================================ *
 * Reads a MCMC config file.					*
 * ============================================================ */
void read_config_mcmc_file(config_mcmc *config, const char *cname, error **err)
{
   FILE *F;
   config_element c = {0, 0.0, ""};
   int i, j;
   char s[500];


   F = fopen_err(cname, "r",err);
   forwardError(*err,__LINE__,);


   /* ============================================================ *
    * Basic sections (Parameters, data)				   *
    * ============================================================ */

   read_config_base(&(config->base), F, 0, err);
   forwardError(*err,__LINE__,);


   /* ============================================================ *
    * MCMC section						   *
    * ============================================================ */


   /* Number of trial points in chain */
   CONFIG_READ(config, nchain, i, F, c, err);

   /* Number of parameter vectors for covariance */
   CONFIG_READ(config, ncov, i, F, c, err);
   if ((config->ncov<=0 || config->ncov>config->nchain) && config->nchain!=0) {
      fprintf(stderr, "Warning: ncov out of range, setting ncov to nchain = %d\n", config->nchain);
      config->ncov = config->nchain;
   }

   /* Burn-in phase = fburnin*ncov */
   CONFIG_READ(config, fburnin, d, F, c, err);

   CONFIG_READ(config, ndecorr, d, F, c, err);

   /* Correction factor f, cov is multiplied with f^2/npar */
   CONFIG_READ(config, fudge, d, F, c, err);

   /* Initial proposal */
   CONFIG_READ_S(&config->base, sinitial, s, F, c, err);
   STRING2ENUM(config->base.initial, config->base.sinitial, initial_t, sinitial_t, j, Ninitial_t, err);
   /* If diagonal, the variance is (max-min)/boxdiv */
   if (config->base.initial==mcmcini_diag) {
      CONFIG_READ(config, boxdiv, d, F, c, err);
   }

   /* Starting point */
   CONFIG_READ_S(config, sstart, s, F, c, err);
   STRING2ENUM(config->start, config->sstart, start_t, sstart_t, j, Nstart_t, err);

   if (config->start==start_fid) {
      config->fid = malloc_err(sizeof(double)*config->base.npar, err);
      forwardError(*err, __LINE__,);
      CONFIG_READ_ARR(config, fid, d, i, config->base.npar, s, F, c, err);
      for (i=0; i<config->base.npar; i++) {
	 testErrorRetVA(config->fid[i]<config->base.min[i], mk_range,
			"Parameter #%d (%g) out of range given by [%g;%g]", *err, __LINE__,,
			i, config->fid[i], config->base.min[i], config->base.max[i]);

	 testErrorRetVA(config->fid[i]>config->base.max[i], mk_range,
			"Parameter #%d (%g) out of range given by [%g;%g]", *err, __LINE__,,
			i, config->fid[i], config->base.min[i], config->base.max[i]);
      }
   } else {
      config->fid = NULL;
   }


   /* Number of bins in histogram */
   CONFIG_READ(&config->base, nbinhist, i, F, c, err);

   fclose(F);
}

void read_config_pmc_file(config_pmc *config, const char *cname, mix_mvdens **proposal, const gsl_rng *rng,
			  int no_init, error **err)
{
   FILE *F;
   config_element c = {0, 0.0, ""};
   int j;


   F = fopen_err(cname, "r", err);
   forwardError(*err, __LINE__,);
 

   /* ============================================================ *
    * Basic sections (Parameters, data)				   *
    * ============================================================ */

   read_config_base(&(config->base), F, no_init, err);
   forwardError(*err,__LINE__,);


   /* ============================================================ *
    * PMC section						   *
    * ============================================================ */

   /* Number of sample points */
   CONFIG_READ(config, nsamples, i, F, c, err);

   /* Multiplicator for final sample */
   CONFIG_READ(config, fsfinal, d, F, c, err);

   /* Number of iterations */
   CONFIG_READ(config, niter, i, F, c, err);

   /* The nclipw largest weights are clipped */
   CONFIG_READ(config, nclipw, i, F, c, err);

   /* === Tempering (new in v1.3) === */
   if (config->base.version >= 1.3) {
      
      CONFIG_READ_S(config, stempering, s, F, c, err);
      STRING2ENUM(config->tempering, config->stempering, tempering_t, stempering_t, j, Ntempering_t, err);
      switch (config->tempering) {
         case tempering_none :
	    break;
         case tempering_linear :
            CONFIG_READ(config, t_min, d, F, c, err);
            break;
         default:
            addErrorVA(mcmc_unknown, "Unknown tempering type (%s)", *err, __LINE__, config->stempering);
            return;
      }

   }

   /* === Proposal === */

   /* Degrees of freedom (-1: mv normal) */
   CONFIG_READ(config, df, i, F, c, err);

   /* Number of components */
   CONFIG_READ(config, ncomp, i, F, c, err);

   /* What to do with dead components */
   CONFIG_READ_S(config, sdead_comp, s, F, c, err);
   STRING2ENUM(config->dead_comp, config->sdead_comp, dead_comp_t, sdead_comp_t, j, Ndead_comp_t, err);

   /* Type of initial proposal */
   CONFIG_READ_S(&config->base, sinitial, s, F, c, err);
   STRING2ENUM(config->base.initial, config->base.sinitial, initial_t, sinitial_t, j, Ninitial_t, err);

   read_initial_proposal_from_config(F, config, proposal, rng, no_init, err); 
   forwardError(*err, __LINE__,);

   /* === Histogram === */

   /* Number of bins in histogram */
   CONFIG_READ(&config->base, nbinhist, i, F, c, err);

   fclose(F);
}

/* ============================================================ *
 * Reads a max_post config file.				*
 * ============================================================ */
void read_config_max_file(config_max *config, const char *cname, error **err)
{
   FILE *F;
   config_element c = {0, 0.0, ""};
   int i, j;
   char s[500];


   F = fopen_err(cname, "r", err);
   forwardError(*err,__LINE__,);


   /* ============================================================ *
    * Basic sections (Parameters, data)				   *
    * ============================================================ */

   read_config_base(&(config->base), F, 0, err);
   forwardError(*err,__LINE__,);


   /* ============================================================ *
    * max_post section						   *
    * ============================================================ */

   /* Starting point */
   CONFIG_READ_S(config, sstart, s, F, c, err);
   STRING2ENUM(config->start, config->sstart, start_t, sstart_t, j, Nstart_t, err);

   if (config->start==start_fid) {
      config->fid = malloc_err(sizeof(double)*config->base.npar, err);
      forwardError(*err, __LINE__,);
      CONFIG_READ_ARR(config, fid, d, i, config->base.npar, s, F, c, err);
      for (i=0; i<config->base.npar; i++) {
	 testErrorRetVA(config->fid[i]<config->base.min[i], mk_range,
			"Parameter #%d (%g) out of range given by [%g;%g]", *err, __LINE__,,
			i, config->fid[i], config->base.min[i], config->base.max[i]);
	 testErrorRetVA(config->fid[i]>config->base.max[i], mk_range,
			"Parameter #%d (%g) out of range given by [%g;%g]", *err, __LINE__,,
			i, config->fid[i], config->base.min[i], config->base.max[i]);
      }
   } else {
      config->fid = NULL;
   }

   CONFIG_READ(config, tolerance, d, F, c, err);
   forwardError(*err, __LINE__,);

   fclose(F);
}

/* ============================================================ *
 * Reads a Fisher matrix config file.				*
 * ============================================================ */
void read_config_fish_file(config_fish *config, const char *cname, error **err)
{
   FILE *F;
   config_element c = {0, 0.0, ""};
   int i;
   char s[500];


   F = fopen_err(cname, "r",err);
   forwardError(*err,__LINE__,);


   /* ============================================================ *
    * Basic sections (Parameters, data)				   *
    * ============================================================ */

   read_config_base(&(config->base), F, 0, err);
   forwardError(*err,__LINE__,);


   /* ============================================================ *
    * Fisher matrix section					   *
    * ============================================================ */

   config->fid = malloc_err(sizeof(double)*config->base.npar, err);
   forwardError(*err, __LINE__,);
   CONFIG_READ_ARR(config, fid, d, i, config->base.npar, s, F, c, err);
   for (i=0; i<config->base.npar; i++) {
      testErrorRetVA(config->fid[i]<config->base.min[i], mk_range,
		     "Parameter #%d (%g) out of range given by [%g;%g]", *err, __LINE__,,
		     i, config->fid[i], config->base.min[i], config->base.max[i]);
      testErrorRetVA(config->fid[i]>config->base.max[i], mk_range,
		     "Parameter #%d (%g) out of range given by [%g;%g]", *err, __LINE__,,
		     i, config->fid[i], config->base.min[i], config->base.max[i]);
   }

   /* Initial proposal */
   CONFIG_READ_S(&config->base, sinitial, s, F, c, err);
   if (strcmp(config->base.sinitial, "Fisher_diag")==0) {
      config->base.initial = mcmcini_fisher_diag;
   } else if (strcmp(config->base.sinitial, "Fisher")==0) {
      config->base.initial = mcmcini_fisher;
   } else {
      *err = addErrorVA(mcmc_unknown, "Wrong Fisher type %s (sinitial in config files), has to be\n"
			"'Fisher' or 'Fisher_diag'", *err, __LINE__, config->base.sinitial);
      return;
   }

}

#define Ntrymax 50
void comp_shift_mean(const config_pmc *config, double *mean, const double *pos0, const gsl_rng *rng,
		     double fshift, error **err)
{
   int i, k;
   double r, boxsize;

   for (i=0,k=0; i<config->base.npar; i++) {
      boxsize = (config->base.max[i]-config->base.min[i]);
      do {
	    r = gsl_ran_flat(rng, -fshift*boxsize, fshift*boxsize);
	    mean[i] = pos0[i] + r;
	    k++;
	    testErrorRet(k>Ntrymax, pmc_tooManySteps,
		        "Could not place proposal component in box. Make sure that fshift "
		        "is not too large and the Fisher matrix is consistent", *err, __LINE__,);
      } while (mean[i]<config->base.min[i] || mean[i]>config->base.max[i]);
   }
}
#undef Ntrymax

void read_initial_proposal_from_config(FILE *F, config_pmc *config, mix_mvdens **proposal_return,
				       const gsl_rng *rng, int no_init, error **err)
{
   config_element c = {0, 0.0, ""};
   int i, j, k, m, sign;
   unsigned long int *sind;
   mvdens *comp=NULL, *fisher;
   mix_mvdens *proposal;
   FILE *FISH, *PROP;
   double r, boxsize, f;
   double *eigenvalue, **eigenvector;
  
   proposal = NULL;

   switch (config->base.initial) {

      case pmcini_none :

	 return;

      case pmcini_fisher_rshift :
      
	 CONFIG_READ(config, fshift, d, F, c, err);
	 CONFIG_READ(config, fvar, d, F, c, err);
      
	 if (no_init==1 || proposal_return==NULL) return;

	 proposal = mix_mvdens_alloc(config->ncomp, config->base.npar, err);     forwardError(*err, __LINE__,);

	 FISH = fopen_err("fisher", "r",err);                   forwardError(*err,__LINE__,);
	 fisher = mvdens_dwnp(FISH, err);                       forwardError(*err, __LINE__,);
	 fclose(FISH);
	 testErrorRetVA(config->base.npar!=fisher->ndim, mk_npar,
			"Fisher matrix dimension (%d) different from config file (%d)",
			*err, __LINE__,, fisher->ndim, config->base.npar);
	 mvdens_inverse(fisher, err);
	 forwardError(*err, __LINE__,);
       
	 /* First component on maximum, will be shifted below */
	 mvdens_from_meanvar(proposal->comp[0], fisher->mean, fisher->std, 1.0);
       
	 comp = mvdens_alloc(config->base.npar, err);            forwardError(*err, __LINE__,);

	 /* Shift and broaden subsequent ncomp components */
	 for (j=0; j<config->ncomp; j++) {

	    /* Mean */
	    comp_shift_mean(config, comp->mean, fisher->mean, rng, config->fshift, err);
        
	    /* Covariance */
	    /* Re-scaling matrix elements */
	    for (i=0; i<config->base.npar; i++) {
	       for (k=0; k<config->base.npar; k++) {
		  comp->std[i*config->base.npar+k] = fisher->std[i*config->base.npar+k]*config->fvar;
	       }
	    }
        
	    mvdens_from_meanvar(proposal->comp[j], comp->mean, comp->std, 1.0);
        
	 }
	 mvdens_free(&fisher);
	 mvdens_free(&comp);
	 break;
      
      case pmcini_file :
      
	 CONFIG_READ_S(config, prop_ini_name, s, F, c, err);
      
	 if (no_init==1 || proposal_return==NULL) return;
      
	 PROP = fopen_err(config->prop_ini_name, "r",err);
	 forwardError(*err,__LINE__,);
	 proposal = mix_mvdens_dwnp(PROP, err);       forwardError(*err, __LINE__,);

	 testErrorRetVA(proposal->ndim!=config->base.npar, pmc_incompat,
			"Initial proposal has incompatible dimension (%d, expected are %d)",
			*err, __LINE__,, proposal->ndim, config->base.npar);
	 testErrorRetVA(proposal->ncomp!=config->ncomp, pmc_incompat,
			"Initial proposal has incompatible number (%d) of components (ncomp = %d in config file)",
			*err, __LINE__,, proposal->ncomp, config->ncomp);
	 for (j=0; j<config->ncomp; j++) {
	    testErrorRetVA(proposal->comp[j]->df!=config->df, pmc_incompat,
			   "Initial proposal has incompatible degrees of freedom %d (df = %d in config file)",
			   *err, __LINE__,, proposal->comp[j]->df, config->df);
	 }
	 break;
      
      case pmcini_fisher_eigen :
      
	 CONFIG_READ(config, fshift, d, F, c, err);
	 CONFIG_READ(config, fvar, d, F, c, err);

	 if (no_init==1 || proposal_return==NULL) return;
      
	 proposal = mix_mvdens_alloc(config->ncomp, config->base.npar, err);     forwardError(*err, __LINE__,);
      
	 FISH = fopen_err("fisher", "r",err);
	 forwardError(*err,__LINE__,);
	 fisher = mvdens_dwnp(FISH, err);                   forwardError(*err, __LINE__,);
	 fclose(FISH);
	 mvdens_inverse(fisher, err);
	 forwardError(*err, __LINE__,);
      
	 eigenvalue  = sm2_vector(1, config->base.npar, err);                    forwardError(*err, __LINE__,);
	 eigenvector = sm2_matrix(1, config->base.npar, 1, config->base.npar, err);        forwardError(*err, __LINE__,);
	 jacobi_transform(fisher->std, config->base.npar, eigenvalue, eigenvector, &m, err); 
	 forwardError(*err, __LINE__,);
      
	 /* Sort eigenvalues */
	 sind  = malloc_err(sizeof(unsigned long int)*(config->base.npar+1), err);   forwardError(*err, __LINE__,);
	 indexx(config->base.npar, eigenvalue, sind, err);				 forwardError(*err, __LINE__,);

	 /* First component on ML point */
	 //mvdens_from_meanvar(proposal->comp[0], fisher->mean, fisher->std, 1.0);

	 comp = mvdens_alloc(config->base.npar, err);				 forwardError(*err, __LINE__,);
	 for (j=0; j<config->base.npar; j++) {
	    for (m=0; m<config->base.npar; m++) {
	       comp->std[j*config->base.npar+m] = config->fvar*fisher->std[j*config->base.npar+m];
	    }
	 }

	 /* Shift components along principal directions */
	 /* i: eigendirection, starting with largest
	    j: dimension
	    k: proposal component
	 */
	 i = config->base.npar;
	 k = 0;
	 while (1) {

	    for (sign=-1; sign<=1; sign+=2) {
	       if (k>=config->ncomp) goto end_shift;
	       for (j=0; j<config->base.npar; j++) {
		  testErrorRet(i<=0, mk_negative, 
			       "Number of proposal components cannot be larger than twice the number of parameters.\n"
			       "For more compoents, set the 'sinitial' flag in the config file to 'fisher_rshift'",
			       *err, __LINE__,);
		  r = config->fshift
		    *eigenvalue[sind[i]]/eigenvalue[sind[config->base.npar]]
		    *eigenvector[sind[i]][j+1];
		  //if (1) fprintf(stderr, "Shift dir=%d comp=%d dim=%d sign=%d   r=%g eval=%g evec=%g\n",
		  //i, k, j, sign, r, eigenvalue[sind[i]], eigenvector[sind[i]][j+1]);
		  comp->mean[j] = fisher->mean[j] + sign*r;
	       }
	       mvdens_from_meanvar(proposal->comp[k], comp->mean, comp->std, 1.0);
	       k++;
	    }
        
	    i--;

	    /* Smallest eigen direction reached? -> Error, too many components */
	    //if (i==0) i = config->base.npar;
	 }
      end_shift:
      
	 free(sind);
	 mvdens_free(&fisher);
	 mvdens_free(&comp);
	 sm2_free_vector(eigenvalue, 1, config->base.npar);
	 sm2_free_matrix(eigenvector, 1, config->base.npar, 1, config->base.npar);

	 break;

      case pmcini_random_pos:

	 CONFIG_READ(config, fmin, d, F, c, err);

	 if (no_init==1 || proposal_return==NULL) return;

	 proposal = mix_mvdens_alloc(config->ncomp, config->base.npar, err);     forwardError(*err, __LINE__,);

	 for (j=0; j<config->ncomp; j++) {

	    for (i=0; i<config->base.npar; i++) {

	       /* Random variance in [fmin;0.5*boxsize] */
	       boxsize = config->base.max[i]-config->base.min[i];
	       r = gsl_ran_flat(rng, config->fmin*boxsize/2.0, 0.5*boxsize/2.0);
	       proposal->comp[j]->std[i*config->base.npar+i] = r*r*r*r;

	       /* Random position in box (near center, depending on variance) */
	       proposal->comp[j]->mean[i] = (config->base.min[i]+config->base.max[i])/2.0;
	       /* The r*FACTOR is ad hoc. The larger the factor, the more central *
		* will the components be placed.					*/
	       f = (boxsize/2.0 - 2*r);
	       if (f>0) proposal->comp[j]->mean[i] += gsl_ran_flat(rng, -1, 1)*f;
	    }
	   
	 }

	 break;

      default :
	 *err = addError(pmc_undef, "Unknown type for initial proposal. Check sinitial (in config file) "
			 "and param.h:sinitial_t", *err, __LINE__);
	 return;
      
   }
  
   if (config->base.initial!=pmcini_file) {
      for (j=0; j<config->ncomp; j++) {
	 proposal->comp[j]->df = config->df;
	 proposal->wght[j] = 1.0/config->ncomp;
      }
   }
  
   //mix_mvdens_print(stderr, proposal);
   mix_mvdens_cholesky_decomp(proposal, err); forwardError(*err, __LINE__,);
  
   *proposal_return = proposal;
}

/* ============================================================ *
 * Reads a data entry from the config file data section and     *
 * Initialises the corresponding likelihood unless no_init=1.   *
 * ============================================================ */
void read_extra_for_data(FILE *F, data_t data, common_like *like, int no_init, error **err)
{
   /* Read data-specific information from config file */
   like->func_wrapper->func_read(&like->state, F, err);
   forwardError(*err, __LINE__,);

   if (data==CMB) {
      if (no_init==0) {
	 check_scamb(like->state, err);
      }
   }

   /* Initialise likelihood structure */
   if (no_init==0) {
      like->func_wrapper->func_init(like, err);
      forwardError(*err, __LINE__,);
   }
}

void out_config_base(FILE *OUT, config_base config, error **err)
{
   int i;

   fprintf(OUT, "out_config_base:\n");
   fprintf(OUT, "===========\n");
   fprintf(OUT, "version    = %g\n", config.version);
   fprintf(OUT, "npar       = %d\n", config.npar);
   fprintf(OUT, "n_ded      = %d\n", config.n_ded);
   fprintf(OUT, "(s)par     = ");
   for (i=0; i<config.npar; i++) fprintf(OUT, " (%s)%d", spar_t(config.par[i]), config.par[i]);
   fprintf(OUT, "\nmin        = ");
   for (i=0; i<config.npar; i++) fprintf(OUT, " % f", config.min[i]);
   fprintf(OUT, "\nmax        = ");
   for (i=0; i<config.npar; i++) fprintf(OUT, " % f", config.max[i]);
   fprintf(OUT, "\nlogpr_def  = %g", config.logpr_default);

   fprintf(OUT, "\nndata      = %d\n", config.ndata);
   fprintf(OUT, "****************\n");
   for (i=0; i<config.ndata; i++) {
      fprintf(OUT, "(s)data[%d]  = (%s)%d\n", i, sdata_t(config.data[i]), config.data[i]);
      if (config.data_extra[i]->func_wrapper->func_print!=NULL) {
	 config.data_extra[i]->func_wrapper->func_print(OUT, config.data_extra[i]->state, err);
	 forwardError(*err, __LINE__,);
      }
   }
   fprintf(OUT, "****************\n");

   fprintf(OUT, "sprior     = %s ", config.sprior);
   if (config.prior!=NULL) {
      fprintf(OUT, "%d", config.nprior);
      for (i=0; i<config.npar; i++) fprintf(OUT, " %d", config.indprior[i]);
   }
   fprintf(OUT, "\n");
}

void out_config_mcmc(FILE *OUT, config_mcmc config, error **err)
{
   int i;

   out_config_base(OUT, config.base, err);
   forwardError(*err, __LINE__,);

   fprintf(OUT, "nchain     = %d\n", config.nchain);
   fprintf(OUT, "ncov       = %d\n", config.ncov);
   fprintf(OUT, "fburnin    = %g\n", config.fburnin);
   fprintf(OUT, "ndecorr    = %d\n", config.ndecorr);
   fprintf(OUT, "fudge      = %g\n", config.fudge);
   fprintf(OUT, "(s)initial = (%s)%d\n", config.base.sinitial, config.base.initial);
   fprintf(OUT, "(s)start   = (%s)%d\n", config.sstart, config.start);
   if (config.start==start_fid) {
      fprintf(OUT, "fid        = ");
      for (i=0; i<config.base.npar; i++) fprintf(OUT, " % f", config.fid[i]);
      fprintf(OUT, "\n");
   }
   fprintf(OUT, "nbinhist   = %d\n", config.base.nbinhist);
}

void out_config_pmc(FILE *OUT, config_pmc config, error **err)
{
   out_config_base(OUT, config.base, err);
   forwardError(*err, __LINE__,);

   fprintf(OUT, "nsamples     = %d\n", config.nsamples);
   fprintf(OUT, "fsfinal      = %g\n", config.fsfinal);
   fprintf(OUT, "niter        = %d\n", config.niter);
   fprintf(OUT, "nclipw       = %d\n", config.nclipw);
   fprintf(OUT, "df           = %d\n", config.df);
   fprintf(OUT, "ncomp        = %d\n", config.ncomp);
   fprintf(OUT, "(s)dead_comp = (%s)%d\n", config.sdead_comp, config.dead_comp);
   fprintf(OUT, "(s)initial   = (%s)%d\n", config.base.sinitial, config.base.initial);

   // TODO: out proposal ...

   fprintf(OUT, "nbinhist     = %d\n", config.base.nbinhist);
}

void out_config_max(FILE *OUT, config_max config, error **err)
{
   int i;

   out_config_base(OUT, config.base, err);
   forwardError(*err, __LINE__,);

   fprintf(OUT, "(s)start   = (%s)%d\n", config.sstart, config.start);
   if (config.start==start_fid) {
      fprintf(OUT, "fid        = ");
      for (i=0; i<config.base.npar; i++) fprintf(OUT, " % f", config.fid[i]);
      fprintf(OUT, "\n");
   }
   fprintf(OUT, "tolerance  = %g\n", config.tolerance);
}

void out_config_fish(FILE *OUT, config_fish config, error **err)
{
   int i;

   out_config_base(OUT, config.base, err);
   forwardError(*err, __LINE__,);

   fprintf(OUT, "fid        = ");
   for (i=0; i<config.base.npar; i++) fprintf(OUT, " % f", config.fid[i]);
   fprintf(OUT, "\n");
   fprintf(OUT, "(s)initial  = (%s)%d\n", config.base.sinitial, config.base.initial);
}

void write_config_base_file(FILE *OUT, config_base config, error **err)
{
   int i;

   fprintf(OUT, "### Basic part ###\n\n");

   fprintf(OUT, "### Parameter section ###\n");
   fprintf(OUT, "npar          %d\n", config.npar);
   fprintf(OUT, "n_ded         %d\n", config.n_ded);
   fprintf(OUT, "spar         ");
   for (i=0; i<config.npar; i++) fprintf(OUT, " %s", spar_t(config.par[i]));
   fprintf(OUT, "\n");
   fprintf(OUT, "min         ");
   for (i=0; i<config.npar; i++) fprintf(OUT, " % g", config.min[i]);
   fprintf(OUT, "\nmax         ");
   for (i=0; i<config.npar; i++) fprintf(OUT, " % g", config.max[i]);

   fprintf(OUT, "\n\n### Data section ###\n");
   fprintf(OUT, "ndata         %d\n", config.ndata);
   for (i=0; i<config.ndata; i++) {
      fprintf(OUT, "sdata         %s\n", sdata_t(config.data[i]));
      forwardError(*err, __LINE__,);
   }
   fprintf(OUT, "\n");

   for (i=0; i<config.ndata; i++) {
      // TODO
      //extra_write_to_config_file(config.data[i], OUT, config.data_extra[i], err);
   }

   fprintf(OUT, "\n### Prior section ###\n");
   fprintf(OUT, "sprior        %s\n", config.sprior);
   if (config.prior!=NULL) {
      fprintf(OUT, "nprior        %d\n", config.nprior);
      fprintf(OUT, "indprior      ");
      for (i=0; i<config.npar; i++)
	fprintf(OUT, "%d", config.indprior[i]);
      fprintf(OUT, "\n");
   }
   fprintf(OUT, "\n");
}

void write_config_pmc_file(FILE *OUT, config_pmc config, error **err)
{
   write_config_base_file(OUT, config.base, err);
   forwardError(*err, __LINE__,);

   fprintf(OUT, "### PMC part ###\n\n");
   fprintf(OUT, "nsamples      %d\n", config.nsamples);
   fprintf(OUT, "fsfinal       %g\n", config.fsfinal);
   fprintf(OUT, "niter         %d\n", config.niter);
   fprintf(OUT, "nclipw        %d\n", config.nclipw);
   fprintf(OUT, "df            %d\n", config.df);
   fprintf(OUT, "ncomp         %d\n", config.ncomp);
   fprintf(OUT, "sdead_comp    %s\n", config.sdead_comp);
   fprintf(OUT, "sinitial      %s\n", config.base.sinitial);

   switch (config.base.initial) {

      case pmcini_fisher_rshift : case pmcini_fisher_eigen :
	 fprintf(OUT, "fshift        %g\n", config.fshift);
	 fprintf(OUT, "fvar          %g\n", config.fvar);
	 break;
      case pmcini_file :
	 fprintf(OUT, "prop_ini_name  %s\n", config.prop_ini_name);
	 break;
      case pmcini_random_pos :
	 fprintf(OUT, "fmin          %g\n", config.fmin);
	 break;
      case pmcini_none :
	 break;
      default :
	 *err = addErrorVA(mk_undef, "Invalid sinitial key '%s(%d)'", *err, __LINE__, 
			   config.base.sinitial, config.base.initial);
	 return;

   }

   fprintf(OUT, "nbinhist      %d\n", config.base.nbinhist);
}

parabox *parabox_from_config(int npar, const double *min, const double *max, error **err)
{
   parabox *pb;
   int i;

   pb = init_parabox(npar, err);   forwardError(*err, __LINE__, NULL);
   for (i=0; i<npar; i++) {
      add_slab(pb, i, min[i], max[i], err);
      forwardError(*err, __LINE__, NULL);
   }
   return pb;
}

/* ============================================================ *
 * The following functions deal with the log-posterior,         *
 * log L = -0.5*chi2 + log(prior).				*
 * Remark: void *extra not const because of mvdens_log_pdf,     *
 * cholesky_decomposition.					*
 * ============================================================ */

double likelihood_log_pdf_single(void *extra, const double *x, error **err)
{
   double res;
   common_like *like;

   like = (common_like*)extra;
 
   res = like->func_wrapper->func_likeli(like, x, err);
   forwardError(*err, __LINE__, MC_LOGBIG);

   return res;
}

double posterior_log_pdf_common_void(void *config, const double *x, error **err)
{
   double res;
   res = posterior_log_pdf_common((config_base*)config, x, err);
   forwardError(*err, __LINE__, -1.0);
   return res;
}



double posterior_log_pdf_common(config_base *config, const double *x, error **err)
{
   double logl, logpost, logpr, pr, sigma_8=-1.0;
   double *xtmp;
   const double *xprior;
   int i, j;
   special_t special;
   common_like *like;


   for (i=0,logpost=0.0; i<config->ndata; i++) {

      /* Before lensing likelihood: get sigma_8 (from CAMB)  *
       * CMB has to come before lensing in the data set list */
      like = (common_like*)(config->data_extra[i]);
      if (config->data[i]==Lensing && sigma_8!=-1) ((lens_state*)(like->state))->sigma_8 = sigma_8;

      logl = likelihood_log_pdf_single(config->data_extra[i], x, err);
      forwardError(*err, __LINE__, 0);
      
      /* After CMB likelihood: save sigma_8 */
      if (config->data[i]==CMB) sigma_8 = ((wmap_state*)(like->state))->sigma_8;

      /* MKDEBUG new v1.4. Prior only added once, not repeatedly for all data sets *
       * Bug fix thanks to R. Schuhmann, B. Joachimi, H. Peiris                    */
      logpost += logl;

      /* Prior. To be added for each data set. */
   }

   /* Prior. To be added once (new v 1.4) */
   if (config->data[0] != Mring) {
      logpr = config->logpr_default;
   } else {
      logpr = 0.0;
   }

   /* Adds special priors for corresponding parameters and removes the default (flat) *
    * prior for those same parameters.                                                */
   // MKDEBUG TODO: define special not for each data set, but globally.
   special = get_special(config->data[0], config->data_extra[0], err);
   forwardError(*err, __LINE__, 0.0);
   logpr += prior_log_pdf_special(special, config->par, config->min, config->max, config->npar, err);
   forwardError(*err, __LINE__, 0.0);

   logpost += logpr;


   /* Additional prior (e.g. previous experiment) */
   // MKDEBUG TODO: Check whether joint likelihood (data * previous exp)
   // is still normalized
   if (config->prior!=NULL) {
      if (config->nprior==0) {
         xprior = x;
      } else {
         xtmp = sm2_vector(0, config->nprior-1, err);
         forwardError(*err, __LINE__, 0);
         for (i=j=0; i<config->npar; i++) {
            if (config->indprior[i]==1) {
               xtmp[j] = x[i];
               j++;
            }
         }
         xprior = (const double*)xtmp;
      }
      pr = mvdens_log_pdf(config->prior, xprior, err);
      forwardError(*err, __LINE__, 0);
      logpost += pr;
   }

   /* Debug info if CMB */
   for (i=0; i<config->ndata; i++) {
      if (config->data[i]==CMB) {
#ifdef COMM_DEBUG
         fprintf(stderr, "wmap logL %g ", logpost); 
         for (j=0; j<config->npar; j++) fprintf(stderr, "%g ", x[j]);
         fprintf(stderr, "\n"); fflush(stderr);
#endif
      }
   }

   return logpost;

}

/* ============================================================ *
 * In case an additional ('special') prior is defined, the      *
 * following is returned:					*
 * log(B_i) - log(S_i)						*
 * where B_i is the default prior ('Box') for the parameter(s)  *
 * i which are affected by the special prior, and S_i is the    *
 * volume of the special prior.					*
 * Since the default prior is -log(B), the result can just be   *
 * added to the default prior to take into account the special  *
 * prior.							*
 * ============================================================ */

double prior_log_pdf_special(special_t special, const par_t *par, const double *min,
      const double *max, int npar, error **err)
{
   double logpr;
   int i, ind_w0, ind_w1;

   switch (special) {

      case none : return 0.0;

      case unity :
	 /* Return - logpr_default */
	 for (i=0,logpr=0.0; i<npar; i++) {
	    logpr += log(max[i] - min[i]);
	 }
	 return logpr;

      case de_conservative :

	 for (i=0,ind_w0=-1,ind_w1=-1; i<npar; i++) {
	    if (par[i]==p_w0de) ind_w0 = i;
	    if (par[i]==p_w1de) ind_w1 = i;
	 }
	 if (ind_w0==-1) return 0.0;  /* One could also print an error message */

	 /* If the w0_de-range is smaller than de_conservative prior, things get more complicated... *
	  * I can't get bothered right now with that.						     */
	 testErrorRetVA(min[ind_w0]>-1 || max[ind_w0]<-1.0/3.0, mk_special_prior,
			"Range of w0_de [%g;%g] should not be smaller than de_conservative prior [-1;-1/3]",
			*err, __LINE__, 0.0, min[ind_w0], max[ind_w0]);

	 if (ind_w1==-1) { /* Only w0_de in parameter list */
	    /* Old volume minus 'special' prior volume */
	    logpr = log(max[ind_w0] - min[ind_w0]) - log(2.0/3.0);
	 } else {
	    /* Rectangle minues two triangle volumes */
	    logpr = log(max[ind_w0] - min[ind_w0]) + log(max[ind_w1] - min[ind_w1])
	      - log(0.5*2.0/3.0*2.0/3.0/(1.0-a_acc)) - log(0.5*2.0/3.0*2.0/3.0);
	 }
	 return logpr;

      default :
	 *err = addErrorVA(mk_undef, "Unknown special prior type %d. Check 'sspecial' keys in config file",
			   *err, __LINE__, special);
	 return 0.0;
   }
}

int retrieve_ded_single(const void *data_extra, double *x_ded, data_t data)
{
   int set_x_ded;
   common_like *like;
   func_deduced_t *fdeduced;

   set_x_ded = 0;

   like = (common_like*)data_extra;

   fdeduced = like->func_wrapper->func_deduced;

   if (fdeduced != NULL) {
      fdeduced(like, x_ded);
   }

   /*
   switch (data) {
      case CMB :
	 //x_ded[0] = ((wmap_state*)like->state)->sigma_8;
	 set_x_ded = 1;
	 break;

      default :
	 break;
   }
   */

   return set_x_ded;
}

void retrieve_ded(const void *extra, double *x_ded, error **err)
{
   int i;
   config_base *config;

   config = (config_base *)extra;

   for (i=0; i<config->ndata; i++) {
      retrieve_ded_single(config->data_extra[i], x_ded, config->data[i]);
   }
}

/* Prints dimensions and limits of a histogram */
void out_histogram(const nd_histogram *hist, FILE *F)
{
   int i;

   if (F) {
      fprintf(F, "histogram:\n");
      fprintf(F, "==========\n");
      fprintf(F, "ndim = %zd, total = %lg, volume = %lg\n", hist->ndim, hist->total, hist->volume);
      for (i=0; i<hist->ndim; i++) {
	 fprintf(F, "%d %zd %lg %lg %lg\n", i, hist->nbins[i], hist->limits[2*i], hist->limits[2*i+1], hist->stps[i]);
      }
      fprintf(F, "\n");
   }
}

/* Prints the content of a histogram, i.e. the log-likelihood data on a grid */
void out_histogram_data(const nd_histogram *hist, FILE *F, error **err)
{
   int idata, ndata, idim, mul, icur;
   double *binc;

   binc = malloc_err(sizeof(double)*hist->ndim, err);  forwardError(*err, __LINE__,);
   for (idim=0,ndata=1; idim<hist->ndim; idim++) {
      ndata *= hist->nbins[idim];
   }

   for (idata=0; idata<ndata; idata++) {
      for (idim=hist->ndim-1,mul=1; idim>=0; idim--) {
	 icur  = (idata/mul)%hist->nbins[idim];
	 binc[idim] = hist->limits[idim*2] + (icur+0.5)*hist->stps[idim];
	 mul  *= hist->nbins[idim];
      }

      for (idim=0; idim<hist->ndim; idim++) {
	 fprintf(F, "% .5f ", binc[idim]);
      }
      fprintf(F, "%.5e\n", hist->data[idata]);
   }
}

/*
void free_config(config_mcmc config, error **err)
{
   int i;

   if (config.indprior!=NULL) free(config.indprior);
   if (config.prior!=NULL) free(config.prior);

   free(config.min); free(config.max); free(config.fid);
   free(config.par);

   for (i=0; i<config.ndata; i++) {
      //mvdens_free((mvdens**)(config.data_extra[i]));
      fprintf(stderr, "Freeing %d\n", i);
      extra_free(config.data[i], config.data_extra+i, err);
      forwardError(*err, __LINE__,);
      fprintf(stderr, "done\n");
   }
   free(config.data);

}
*/

/* ============================================================ *
 * Creates 1d+2d histograms of all parameter combinations.	*
 * (Calls following two functions.)				*
 * ============================================================ */

void create_1d_2d_hist(double *accstates, long naccepted, int nbinhist, int npar, const double *min,
		    const double *max, double *weights, FILE *FLOG, const char *path, error **err)
{
   create_1d_hist(accstates, naccepted, nbinhist, npar, min, max, weights, path, err);
   checkPrintErr_and_continue(err, stderr, FLOG, NULL, 0);

   create_2d_hist(accstates, naccepted, nbinhist, npar, min, max, weights, path, err);
   checkPrintErr_and_continue(err, stderr, FLOG, NULL, 0);
}

/* ============================================================ *
 * Creates 1d histograms of all parameters.			*
 * ============================================================ */

void create_1d_hist(double *accstates, long naccepted, int nbinhist, int npar, const double *min,
		    const double *max, double *weights, const char *path, sm2_error **err)
{
  nd_histogram *hist;
  size_t *nbins, *pidx;
  double *limits, *w;
  int ndim, p, i;
  FILE *F;
  char name[100];
  
  if (nbinhist==0) return;
  
  ndim     = 1;
  nbins    = malloc_err(sizeof(size_t)*ndim, err);   forwardError(*err, __LINE__,);
  nbins[0] = nbinhist;
  limits   = malloc_err(sizeof(double)*ndim*2, err); forwardError(*err, __LINE__,);
  pidx     = malloc_err(sizeof(size_t)*ndim, err);   forwardError(*err, __LINE__,);

  /* TODO: This is only necessary to get the shading right in likeli.i:plot_chi2n */
  if (weights!=NULL) {
     w = malloc_err(sizeof(double)*naccepted, err); forwardError(*err, __LINE__,);
     for (i=0; i<naccepted; i++) w[i] = weights[i]*naccepted;
  } else{
     w = NULL;
  }

  for (p=0; p<npar; p++) {
    
    limits[0] = min[p];
    limits[1] = max[p];
    
    pidx[0] = p;
    
    if (nbinhist>0) {
       hist = init_nd_histogram(ndim, nbins, limits, err);
       forwardError(*err, __LINE__,);
       acc_histogram(npar, naccepted, accstates, w, pidx, hist, err);
       forwardError(*err, __LINE__,);
    } else {
      hist = make_histogram(ndim, pidx, npar, naccepted, accstates, w, err);
      forwardError(*err, __LINE__,);
    }
    
    //out_histogram(hist, FLOG);
    
    sprintf(name, "./%s/%s_%d", path, hist_name, p);
    F = fopen_err(name, "w",err);
    forwardError(*err,__LINE__,);

    out_histogram_data(hist, F, err);       forwardError(*err, __LINE__,);
    fclose(F);

    free_nd_histogram(&hist);

  }

  free(nbins);
  free(limits);
  free(pidx);
}

/* ============================================================ *
 * Creates 2d histograms of all parameter combinations.		*
 * ============================================================ */

void create_2d_hist(double *accstates, long naccepted, int nbinhist, int npar, const double *min,
		    const double *max, double *weights, const char *path, sm2_error **err)
{
  nd_histogram *hist;
  size_t *nbins, *pidx;
  double *limits, *w;
  int ndim, p, q, i;
  FILE *F;
  char name[100];


  if (nbinhist==0) return;
  
  ndim     = 2;
  nbins    = malloc_err(sizeof(size_t)*ndim, err);   forwardError(*err, __LINE__,);
  nbins[0] = nbins[1] = nbinhist;
  limits   = malloc_err(sizeof(double)*ndim*2, err); forwardError(*err, __LINE__,);
  pidx     = malloc_err(sizeof(size_t)*ndim, err);   forwardError(*err, __LINE__,);
  
  /* TODO: This is only necessary to get the shading right in likeli.i:plot_chi2n */
  if (weights!=NULL) {
     w = malloc_err(sizeof(double)*naccepted, err); forwardError(*err, __LINE__,);
     for (i=0; i<naccepted; i++) w[i] = weights[i]*naccepted;
  } else{
     w = NULL;
  }

  for (p=0; p<npar-1; p++) {
    for (q=p+1; q<npar; q++) {
      
      limits[0] = min[p];
      limits[1] = max[p];
      limits[2] = min[q];
      limits[3] = max[q];
      
      pidx[0] = p;
      pidx[1] = q;
      
      if (nbinhist>0) {
        hist = init_nd_histogram(ndim, nbins, limits, err);
        forwardError(*err, __LINE__,);
        acc_histogram(npar, naccepted, accstates, w, pidx, hist, err);
        forwardError(*err, __LINE__,);
      } else {
        hist = make_histogram(ndim, pidx, npar, naccepted, accstates, w, err);
      }
      
      //out_histogram(hist, FLOG);
      
      sprintf(name, "./%s/%s_%d_%d", path, hist_name, p, q);
      F = fopen_err(name, "w",err);
      forwardError(*err,__LINE__,)
      out_histogram_data(hist, F, err);
      forwardError(*err, __LINE__,);
      fclose(F);
      
      free_nd_histogram(&hist);
    }
  }
  
  free(nbins);
  free(limits);
  free(pidx);
}

void check_npar(int npar, int nmax, error **err)
{
   testErrorRet(nmax>npar, mk_npar, "Number of parameters required for this parametrization too large "
		"(check npar and sparam in config file", *err, __LINE__,);
}

double *merge_param_param_ded(double *param, const double *param_ded, long nsample, int npar, int n_ded, error **err)
{
   int i, d;
   double *Xall;

   if (n_ded>0) {
      Xall = malloc_err(sizeof(double)*nsample*(npar+n_ded), err);
      forwardError(*err, __LINE__, NULL);

      for (i=0; i<nsample; i++) {
	 for (d=0; d<npar; d++) {
	    Xall[i*(npar+n_ded)+d] = param[i*npar+d];
	 }
	 for (d=0; d<n_ded; d++) {
	    Xall[i*(npar+n_ded)+d+npar] = param_ded[i*n_ded+d];
	 }
      }
   } else {
      Xall = param;
   }

   return Xall;
}

void do_nothing()
{
}

/* ============================================================ *
 * The following functions were moved from wrapper.c		*
 * ============================================================ */

common_like *extra_init_common(data_t i, const par_t *par, int npar, error **err)
{
   common_like *like;
   init_functions_t *init_func;

   like = malloc_err(sizeof(common_like), err);
   forwardError(*err, __LINE__, NULL);

   /* TODO: sizeof(void) gives warning */
   like->state = malloc_err(sizeof(void), err);
   forwardError(*err, __LINE__, NULL);

   /* Copy parameter information */
   like->npar = npar;
   like->par  = copy_par_t(par, npar, err);
   forwardError(*err, __LINE__, NULL);

   /* Initialise wrapper functions */
   init_func          = init_func_t(i, err);     forwardError(*err, __LINE__, NULL);
   like->func_wrapper = init_func(err);          forwardError(*err, __LINE__, NULL);

   like->want_model = 0;
   like->model      = NULL;
   like->Nmodel     = 0;

   return like;
}

special_t get_special(data_t data, common_like *like, error **err)
{
   special_t special;

   if (like->func_wrapper->func_special==NULL) {
      return none;
   } else {
      special = like->func_wrapper->func_special(like->state);
   }

   return special;
}

/* ============================================================ *
 * 'Wrapper' functions for mvdens and mix_mvdens.		*
 * ============================================================ */

functions_wrapper_t *init_functions_Mvdens(error **err)
{
   functions_wrapper_t *init_func;

   init_func = init_functions_wrapper(read_from_config_Mvdens, init_dummy, likeli_Mvdens,
				      NULL, NULL, NULL, err);
   forwardError(*err, __LINE__, NULL);

   return init_func;
}

void read_mvdens_component(mvdens *g, FILE *F, error **err)
{
   config_element c = {0, 0.0, ""};
   int i;
   char s[1024];

   CONFIG_READ_ARR(g, mean, d, i, g->ndim, s, F, c, err);
   CONFIG_READ_ARR(g, std, d, i, g->ndim*g->ndim, s, F, c, err);
   CONFIG_READ(g, df, i, F, c, err);
}

void read_from_config_Mvdens(void **state, FILE *F, error **err)
{
   mvdens *g, *tmp;
   config_element c = {0, 0.0, ""};

   tmp = mvdens_alloc(1, err);                      forwardError(*err, __LINE__,);
   CONFIG_READ(tmp, ndim, i, F, c, err);
   g = mvdens_alloc(tmp->ndim, err);                forwardError(*err, __LINE__,);
   mvdens_free(&tmp);

   read_mvdens_component(g, F, err);                forwardError(*err, __LINE__,);

   mvdens_cholesky_decomp(g, err);                  forwardError(*err, __LINE__,);
   *state = g;
}

void init_dummy(common_like *like, error **err)
{
   /* Do nothing */
}

double likeli_Mvdens(common_like *like, const double *params, error **err)
{
   double res;

   res = mvdens_log_pdf_void((void *)(like->state), params, err);
   forwardError(*err, __LINE__, -1.0);

   return res;
}

functions_wrapper_t *init_functions_MixMvdens(error **err)
{
   functions_wrapper_t *init_func;

   init_func = init_functions_wrapper(read_from_config_MixMvdens, init_dummy, likeli_MixMvdens,
				      NULL, NULL, NULL, err);
   forwardError(*err, __LINE__, NULL);

   return init_func;
}

void read_from_config_MixMvdens(void **state, FILE *F, error **err)
{
   mix_mvdens *g, *tmp;
   int k;
   config_element c = {0, 0.0, ""};

   tmp = mix_mvdens_alloc(1, 1, err);               forwardError(*err, __LINE__,);
   CONFIG_READ(tmp, ncomp, i, F, c, err);
   CONFIG_READ(tmp, ndim, i, F, c, err);
   tmp = mix_mvdens_alloc(tmp->ncomp, tmp->ndim, err);
   forwardError(*err, __LINE__,);
   g = mix_mvdens_alloc(tmp->ncomp, tmp->ndim, err); forwardError(*err, __LINE__,);
   mix_mvdens_free(&tmp);

   for (k=0; k<g->ncomp; k++) {
      read_mvdens_component(g->comp[k], F, err);    forwardError(*err, __LINE__,);
      forwardError(*err, __LINE__,);
      g->wght[k] = 1.0/g->ncomp;
   }

   mix_mvdens_cholesky_decomp(g, err);              forwardError(*err, __LINE__,);

   *state = g;
}

double likeli_MixMvdens(common_like *like, const double *params, error **err)
{
   double res;
   res = mix_mvdens_log_pdf_void((void *)(like->state), params, err);
   forwardError(*err, __LINE__, -1.0);
   return res;
}


/* ============================================================ *
 * Sets the four basic density parameters Omegam, Omegab,       *
 * Omegade, Omeganumass in the structure cosmo.			*
 * ============================================================ */
void set_base_parameters(cosmo *cosmo, double Omegam, double Omegab, double Omegade, double Omeganumass,
			 double Omegac, double OmegaK,
			 double omegam, double omegab, double omegade, double omeganumass,
			 double omegac, double omegaK, double h100,
			 int iOmegade, int iOmegaK, int iomegade, int iomegaK, error **err)
{
   if (Omegam>0 || Omegab>0 || iOmegade==1 || Omeganumass>0 || Omegac>0 || iOmegaK==1) {

      /* Check whether all parameters are consistently either physical (omegaX) or non-physical (OmegaX) */
      testErrorRet(omegam>0 || omegab>0 || iomegade==1 || omeganumass>0 || omegac>0 || iomegaK==1, tls_cosmo_par,
		   "Mixing of physical and non-physical density parameters not possible", *err, __LINE__,);

   } else {

      testErrorRetVA(h100<0, tls_cosmo_par, "Hubble parameter negative", *err, __LINE__,, h100);

      /* For the Coyote 2010 emulator, set h according to the CMB peak constraint. */
      if (cosmo->nonlinear == coyote10) {
         /* TODO: make sure, w0_de is already set correctly */
         h100 = getH0fromCMB(omegam, omegab, cosmo->w0_de, 1);
      }

      /* Transform all parameters to non-physical ones */
      Omegam  = omegam/h100/h100;
      Omegab  = omegab/h100/h100;
      Omegade = omegade/h100/h100;
      Omeganumass = omeganumass/h100/h100;
      Omegac  = omegac/h100/h100;
      OmegaK  = omegaK/h100/h100;
      iOmegade = iomegade;
      iOmegaK  = iomegaK;

   }


   /* Check total density consistencies */
   testErrorRet(Omegam>0 && iOmegade==1 && iOmegaK==1, tls_cosmo_par,
		"Total density overdetermined", *err, __LINE__,);
   testErrorRet(Omegam>0 && Omegab>0 && Omegac>0, tls_cosmo_par,
		"Matter density overdetermined", *err, __LINE__,);

   /* Set massive neutrino density to default value zero */
   if (Omeganumass<0) Omeganumass = 0;

   /* Set the 5 base parameters a partir de tous les autres */
   set_base_Omegam(cosmo, Omegam, Omegab, Omegade, Omeganumass, Omegac, OmegaK, iOmegaK);
   set_base_Omegab(cosmo, Omegam, Omegab, Omegac);
   set_base_Omegade(cosmo, Omegam, Omegab, Omegade, Omeganumass, Omegac, OmegaK);
   set_base_Omeganumass(cosmo, Omeganumass);
}

void set_base_Omegam(cosmo *cosmo, double Omegam, double Omegab, double Omegade, double Omeganumass,
		     double Omegac, double OmegaK, int iOmegaK)
{
   if (Omegam>0) {

      cosmo->Omega_m = Omegam;

   } else if (Omegab>0 && Omegac>0) {

      cosmo->Omega_m = Omegab + Omegac;

   } else if (Omegade>0 && iOmegaK==1) {

      cosmo->Omega_m = 1.0 - Omegade - OmegaK - Omeganumass;

   } else {

      /* cosmo->Omega_m has default value */

   }
}

void set_base_Omegab(cosmo *cosmo, double Omegam, double Omegab, double Omegac)
{
   if (Omegab>0) {

      cosmo->Omega_b = Omegab;

   } else if (Omegam>0 && Omegac>0) {

      cosmo->Omega_b = Omegam - Omegac;

   } else {

      /* cosmo->Omega_b has default value */

   }
}

void set_base_Omegade(cosmo *cosmo, double Omegam, double Omegab, double Omegade, double Omeganumass,
		      double Omegac, double OmegaK)
{
   if (Omegade>0) {

      cosmo->Omega_de = Omegade;

   } else if (Omegam>0) {

      /* Assures flat Universe if OmegaK=0 (default) */
      cosmo->Omega_de = 1.0 - Omegam - OmegaK - Omeganumass;

   } else if (Omegab>0 && Omegac>0) {

      /* Assures flat Universe if OmegaK=0 (default) */
      cosmo->Omega_de = 1.0 - Omegab - Omegac - OmegaK - Omeganumass;

   } else {

      /* cosmo->Omega_de has default value */

   }
}

void set_base_Omeganumass(cosmo *cosmo, double Omeganumass)
{
   if (Omeganumass>0) cosmo->Omega_nu_mass = Omeganumass;
}

void reset_base_parameters(double *Omegam, double *Omegab, double *Omegade, double *Omeganumass,
			   double *Omegac, double *OmegaK,
			   double *omegam, double *omegab, double *omegade, double *omeganumass,
			   double *omegac, double *omegaK, double *h100,
			   int *iOmegade, int *iOmegaK, int *iomegade, int *iomegaK)
{
   *Omegam = *Omegab = *Omegac = *Omegade = *Omeganumass = -1.0;
   *omegam = *omegab = *omegac = *omegade = *omeganumass = *h100 = -1.0;
   *OmegaK = *omegaK = 0.0;
   *iOmegade = *iOmegaK = *iomegade = *iomegaK = 0;
}


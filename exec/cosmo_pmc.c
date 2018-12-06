/* ============================================================ *
 * cosmo_pmc.c							*
 * Martin Kilbinger 2008-2010					*
 * ============================================================ */


#include "cosmo_pmc.h"


void clean_previous_run(int iter, int niter, const char iterdirname[])
{
   int i;
   char name[NSTR];

   /* Delete subsequent simulations and proposals (from previous run) */
   for (i=iter+1; i<niter; i++) {
      sprintf(name, "%s_%d/%s", iter_dir, i, pmcsim_name);
      unlink(name);
      sprintf(name, "%s_%d/%s", iter_dir, i, proposal_name);
      unlink(name);
   }
}

void write_proposal(mix_mvdens *proposal, FILE *F1, FILE *F2, int n, const char *iterdirname, error **err)
{
   char cname[NSTR];
   FILE *F;

   if (F1!=NULL) mix_mvdens_print(F1, proposal);
   if (F2!=NULL) mix_mvdens_print(F2, proposal);
   if (n>=0) {
      sprintf(cname, "%s/%s", iterdirname, proposal_name);
   } else {
      sprintf(cname, "%s_%s", proposal_name, "fin");
   }
   F = fopen_err(cname, "w", err);
   forwardError(*err,__LINE__,);
   mix_mvdens_dump(F, proposal);
   fclose(F);
}

void write_perplexity_and_ess(pmc_simu *psim, int iter, int sum_nsamples, double *sum_ess, FILE *PERP, error **err)
{
   double perp, ess;

   perp = perplexity_and_ess(psim, MC_UNORM, &ess, err);
   forwardError(*err, __LINE__,);
   *sum_ess += ess;

   /* Header */
   if (iter==0) fprintf(PERP, "# iter nsample perplexity   ess/it      ess\n");

   /* Writing sum_nsamples ignores clipped points... */
   fprintf(PERP, "    %2d  %6d  %#9.5g %#8.1F %#8.1F\n", iter, sum_nsamples, perp, ess, *sum_ess);
   fflush(PERP);
}

void write_evidence(pmc_simu *psim, int iter, FILE *EVI, error **err)
{
   double evi, ln_evi;

   evi = evidence(psim, &ln_evi, err);
   quitOnError(*err, __LINE__, stderr);

   /* Header */
   if (iter==0) fprintf(EVI, "# iter   log10(E)    ln(E)    E\n");

   fprintf(EVI, "% 6d % g % g % g\n", iter, ln_evi*M_LOG10E, ln_evi, evi);
   fflush(EVI);
}

void write_enc(mix_mvdens *proposal, int iter, FILE *ENC, error **err)
{
   double enc;

   /* Header */
   if (iter==0) fprintf(ENC, "# iter effective_number_of_components\n");

   enc = effective_number_of_components(proposal, err); quitOnError(*err, __LINE__, stderr);
   fprintf(ENC, "    %2d  %g\n", iter, enc);
   fflush(ENC);
}

void histograms_and_covariance(pmc_simu *psim, const char *iterdirname, config_base *config, FILE *FLOG,
			       error **err)
{
   double *Xall;
   char command[1024], cname[NSTR], pname[NSTR];

   /* Remove old files in subdirectory iter_[iter] */
   sprintf(command, "rm -f %s/%s_*", iterdirname, hist_name);
   system(command);

   Xall = merge_param_param_ded(psim->X, psim->X_ded, psim->nsamples, config->npar,
					 config->n_ded, err);
   forwardError(*err, __LINE__,);

   create_1d_2d_hist(Xall, psim->nsamples, config->nbinhist, config->npar+config->n_ded,
		     config->min, config->max, psim->weights, FLOG, iterdirname, err);
   forwardError(*err, __LINE__,);

   covariance_from_sample_all(psim->X, psim->weights, psim->ndim, psim->nsamples,
			      Xall, psim->ndim+config->n_ded, iterdirname, err);
   checkPrintErr_and_continue(err, stderr, FLOG, NULL, 0);
   sprintf(pname, "%s/%s%s.%s", iterdirname, covar_name, inv_suf, chainfin_suf);
   sprintf(cname, "%s/%s_%s%s", iterdirname, evidence_name, covar_name, inv_suf);
   /* Laplace approximation using the inverse covariance */
   evidence_approx(config, pname, cname, err);
   checkPrintErr_and_continue(err, stderr, FLOG, NULL, 0);

   if (config->n_ded>0) free(Xall);
}

/* ============================================================ *
 * Calculates the Laplace approximation of the evidence,        *
 * reading the file 'covname' as inverse covariance or Fisher   *
 * matrix.							*
 * ============================================================ */
void evidence_approx(config_base *config, const char *covname, const char *outname, error **err)
{
   FILE *MVD, *EVI;
   mvdens *mvd;
   double det, logmaxP, ln_evi, evi, volume;
   int i;


   MVD = fopen(covname, "r");
   if (MVD==NULL) return;

   mvd = mvdens_dwnp(MVD, err);                      forwardError(*err, __LINE__,);
   fclose(MVD);

   testErrorRetVA(config->npar!=mvd->ndim, mv_dimension, "Wrong dimension %d in Fisher matrix, expected %d",
		  *err, __LINE__,, mvd->ndim, config->npar);

   //det = sm2_inverse(mvd->std, mvd->ndim, err);   forwardError(*err, __LINE__,);
   det = mvdens_inverse(mvd, err); forwardError(*err, __LINE__,);
   testErrorRet(det<=0, mcmc_negative, "Covariance matrix is not positive definite", *err, __LINE__,);

   logmaxP = posterior_log_pdf_common(config, mvd->mean, err);
   forwardError(*err, __LINE__,);

   ln_evi = 0.5*log(2.0*pi)*config->npar - 0.5*log(det) + logmaxP;
   evi = exp(ln_evi);

   EVI = fopen_err(outname, "w",err); forwardError(*err,__LINE__,);
   //fprintf(EVI, "# log10(E) ln(E) E [Laplace approximation]\n");
   fprintf(EVI, "# iter   log10(E)    ln(E)    E\n");
   fprintf(EVI, "%6d % g % g % g\n", -1, ln_evi*M_LOG10E, ln_evi, evi);

   /* The following is unused, volume factor (=1/prior) is already in posterior */
   for (i=0,volume=0.0; i<config->npar; i++) {
      volume += log(config->max[i] - config->min[i]);
   }
   fprintf(EVI, "# 0.5*ln2pi*n=%g 0.5*log|F|=%g logmaxP=%g (logmaxL=%g volume=%g)\n", 
	   0.5*log(2.0*pi)*config->npar, -0.5*log(det), logmaxP, logmaxP-volume, volume);

   fclose(EVI);
}

#define FSHIFT_DEFAULT 0.1
void revive_comp(const config_pmc *config, mix_mvdens *proposal, int i, int Nrevive, const gsl_rng *rng, error **err)
{
   int j, imax;
   double wmax, fshift;

   /* Look for component with maximum weight */
   for (j=0,imax=-1,wmax=0.0; j<proposal->ncomp; j++) {
      if (proposal->wght[j]>wmax) {
	 wmax = proposal->wght[j];
	 imax = j;
      }
   }
   testErrorRet(imax==-1, math_negative, "Components weights are all zero", *err, __LINE__,);

   if (config->base.initial==pmcini_fisher_rshift || config->base.initial==pmcini_fisher_eigen) {
      fshift = config->fshift;
   } else {
      fshift = FSHIFT_DEFAULT;
   }
   /* Decrease shift parameter if in previous iterations components have been revived already */
   fshift /= (double)Nrevive;

   /* Place new dead component near maximum-weight component */
   mvdens_from_meanvar(proposal->comp[i], proposal->comp[imax]->mean, proposal->comp[imax]->std, 1.0);
   comp_shift_mean(config, proposal->comp[i]->mean, proposal->comp[imax]->mean, rng, fshift, err);
   forwardError(*err, __LINE__,);

   /* Maximum and new component get half of old (maximum) weight */
   proposal->wght[imax] /= 2.0;
   proposal->wght[i]     = proposal->wght[imax];
   testErrorRetVA(!isfinite(proposal->wght[i]), pmc_negWeight, "Proposal weight #%d is %g", *err, __LINE__,, i, proposal->wght[i]);
}

#define EPS 1.0e-12
void update_and_deal_with_dead(const config_pmc *config, mix_mvdens *proposal, pmc_simu *psim,
			       int *Nrevive, const gsl_rng *rng, error **err)
{
   int i;
   double wght_sum;


   /* Rao-Blackwell update */
   update_prop_rb(proposal, psim, err);
   quitOnError(*err, __LINE__, stderr);


   /* Deal with dead components */
   switch (config->dead_comp) {
      case bury :
	 /* Components already declared dead and buried in cleanup_after_update */
	 break;

      case revive :
	 /* Double revive factor to decrease shift for mean from maximum component */
	 *Nrevive *= 2;

	 for (i=0;i<proposal->ncomp;i++) {
	    if (proposal->wght[i]<EPS) {
	       revive_comp(config, proposal, i, *Nrevive, rng, err);
	       forwardError(*err, __LINE__,);
	    }
	 }

	 /* Renormalize weights */
	 wght_sum = 0.0;
	 for (i=0; i<proposal->ncomp; i++) {
	    wght_sum += proposal->wght[i];
	 }
	 testErrorRetVA(wght_sum<=0.0, pmc_negWeight, "Sum of weight is %g", *err, __LINE__,, wght_sum);
	 testErrorRetVA(!isfinite(wght_sum), pmc_negWeight, "Sum of weight is %g", *err, __LINE__,, wght_sum);
	 gsl_vector_scale(proposal->wght_view, 1.0/wght_sum);
	 break;

      default :
	 *err = addErrorVA(mk_undef, "Dead component mode (%d) undefined, set 'sdead_comp' in config file",
			   *err, __LINE__, config->dead_comp);
	 return;
   }

}
#undef EPS


/* ============================================================ *
 * If the pmc sample ('pmcsim') and  proposal files do not      *
 * exist, a PMC iteration is run. In the contrary case, the     *
 * files are read.						                            *
 * ============================================================ */
void run_pmc_iteration_MPI(pmc_simu *psim, mix_mvdens **proposal_p, int iter, int this_nsamples, int *Nrevive, config_pmc *config,
      const char *iterdirname, int myid, int nproc, int quiet, gsl_rng *rng, parabox *pb, FILE *FLOG, error **err)
{
   char cname[NSTR];
   int nok, mysamples, master_samples=0;
   mix_mvdens *proposal, *proposal_local = NULL;
   double norm;

   proposal = *proposal_p;

   /* Sample from proposal, send and receive simulation */
   if (myid==0) {     /* Server */

      fprintf(FLOG, "=== Iteration #%d ===\n", iter);
      if (! quiet) fprintf(stderr, "=== Iteration #%d ===\n", iter);

      clean_previous_run(iter, config->niter, iterdirname);

      /* Increase size in case of final iteration */
      pmc_simu_realloc(psim, this_nsamples, err);
      forwardError(*err, __LINE__,);

      write_proposal(proposal, NULL, FLOG, iter, iterdirname, err);            forwardError(*err, __LINE__,);

      /* Sample from proposal (simulate PMC sample) */
      fprintf(FLOG, "proc %d >> simulating (iter %d)\n", myid, iter);
      nok = simulate_mix_mvdens(psim, proposal, rng, pb, err);                 forwardError(*err, __LINE__,);
      if (nok==0) {
         *err = addError(pmc_nosamplep, "No point simulated", *err, __LINE__);
         return;
      }
      fprintf(FLOG, "proc %d >> sendinf pmc sim (iter %d)\n", myid, iter);
      mysamples = send_simulation(psim, nproc, err);                           forwardError(*err, __LINE__,);
      master_samples = psim->nsamples;
      psim->nsamples = mysamples;

   } else {    /* Clients */

      /* Receive simulation */
      fprintf(FLOG, "proc %d >> receive pmc sim (iter %d)\n", myid, iter);
      receive_simulation(psim, nproc, myid, err);                              forwardError(*err, __LINE__,);

   }

   /* Calculation of the PMC weights (server and clients) */

   fprintf(FLOG, "proc %d >> work on %ld samples (iter %d)\n", myid, psim->nsamples, iter);
   /* Compute importance weights */
   nok = generic_get_importance_weight_and_deduced_verb(psim, proposal, mix_mvdens_log_pdf_void,
         posterior_log_pdf_common_void, retrieve_ded,
         (void*)&(config->base), quiet, err);
   forwardError(*err, __LINE__,);
   fprintf(FLOG, "proc %d >> finished importance weights (iter %d), nok=%d\n", myid, iter, nok);


   /* Communicate importance weights, update and communicate proposal */
   if (myid!=0) {    /* Clients */

      /* Send importance weights */
      fprintf(FLOG, "proc %d >> send weights (iter %d)\n" ,myid, iter);
      send_importance_weight(myid, nproc, psim, nok, err);    forwardError(*err, __LINE__,);

      /* Receive updated proposal */

      fprintf(FLOG, "proc %d >> receive proposal (iter %d)\n", myid, iter);
      proposal_local = receive_mix_mvdens(myid, nproc, err);        forwardError(*err, __LINE__,);

      mix_mvdens_copy(proposal, proposal_local, err);               forwardError(*err, __LINE__,);

      /* The following line gives glibc error/segmentation fault on magique. *
       * As long as this bug (?) is not fixed, the memory is unfortunately   *
       * freed.								     */
      //mix_mvdens_free(&proposal_local);

   } else {     /* Server */

      mysamples      = psim->nsamples;
      psim->nsamples = master_samples;      
      fprintf(FLOG, "proc %d >> receive weights (iter %d)\n", myid, iter);
      nok = receive_importance_weight(psim, nproc, nok, mysamples, err);
      forwardError(*err, __LINE__,);

      fprintf(FLOG, "Server: normalize_importance_weight (iter %d)\n", iter);
      norm = normalize_importance_weight(psim, err);
      forwardError(*err, __LINE__,);

      fprintf(FLOG, "Server: update_prop_rb (iter %d)\n", iter);
      update_and_deal_with_dead(config, proposal, psim, Nrevive, rng, err);
      fprintf(FLOG, "Nrevive after update: %d\n", *Nrevive);

      /* Broadcast results to clients */
      fprintf(FLOG, "proc %d >> send proposal (iter %d)\n", myid, iter);
      send_mix_mvdens(proposal, nproc, err);
      forwardError(*err, __LINE__,);

      /* Write simulation (PMC points, same format as MCM chain) */
      sprintf(cname, "%s/%s", iterdirname, pmcsim_name);
      out_pmc_simu_cosmo_pmc(cname, psim, config->base.par, norm, err);
      forwardError(*err, __LINE__,);

      /* NEW: clip weights here as well. The other call is in pmc_simu_from_file */
      if (config->nclipw>0) {
         clip_weights(psim, config->nclipw, stderr, err);
         forwardError(*err, __LINE__,);
      }

   }
}

pmc_simu *read_pmc_iteration_MPI(mix_mvdens **proposal, int iter, int this_nsamples, config_pmc *config,
      int myid, int nproc, FILE *PMCSIM, FILE *PROP, FILE *FLOG, error **err)
{
   pmc_simu *psim = NULL;

   /* Read PMC simulation and proposal */
   if (myid==0) {    /* server */

      if (*proposal) mix_mvdens_free(proposal);
      *proposal = mix_mvdens_dwnp(PROP, err);
      forwardError(*err, __LINE__, NULL);

      psim = pmc_simu_from_file(PMCSIM, this_nsamples, config->base.npar, config->base.n_ded, *proposal,
            config->nclipw, err);
      forwardError(*err, __LINE__, NULL);

      /* Broadcast results to clients */
      fprintf(FLOG, "proc %d >> send proposal (iter %d)\n", myid, iter);
      send_mix_mvdens(*proposal, nproc, err);
      forwardError(*err, __LINE__, NULL);

   } else {    /* clients */

      /* Init psim container with dummy size */
      psim = pmc_simu_init_plus_ded(1, config->base.npar, config->base.n_ded, err);
      forwardError(*err, __LINE__, NULL);

      /* Receive updated proposal */
      mix_mvdens_free(proposal);                              forwardError(*err, __LINE__, NULL);
      fprintf(FLOG, "proc %d >> receive proposal (iter %d)\n", myid, iter);
      *proposal = receive_mix_mvdens(myid, nproc, err);        forwardError(*err, __LINE__, NULL);

   }

   return psim;
}

void post_processing(pmc_simu *psim, mix_mvdens *proposal, int iter, int sum_nsamples, double *sum_ess, 
      config_base *config, const char *iterdirname, FILE *PERP, FILE *EVI, FILE *ENC,
      FILE *FLOG, error **err)
{
   char pname[NSTR];
   double *pmean;

   write_perplexity_and_ess(psim, iter, sum_nsamples, sum_ess, PERP, err);  forwardError(*err, __LINE__,);
   write_evidence(psim, iter, EVI, err);                    	             forwardError(*err, __LINE__,);
   /* Proposal is actually used in the *next* iteration */
   write_enc(proposal, iter+1, ENC, err);				    forwardError(*err, __LINE__,);

   /* Mean and 1d confidence levels */
   sprintf(pname, "%s/%s", iterdirname, mean_name);
   pmean = mean_sigma_from_psim(psim, pname, config->par, avg_mean, err);    forwardError(*err, __LINE__,);
   sprintf(pname, "%s/%s", iterdirname, model_mean_name);
   write_mean(pmean, pname, config->ndata, (void**)config->data_extra, err); forwardError(*err, __LINE__,);

   /* Histograms (if config->nbinhist!=0) */
   histograms_and_covariance(psim, iterdirname, config, FLOG, err);         forwardError(*err, __LINE__,);
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

   fprintf(stderr, "Usage: cosmo_pmc [OPTIONS]\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: 'config_pmc')\n");
   fprintf(stderr, "  -s SEED          Use SEED for random number generator. If SEED=-1 (default)\n");
   fprintf(stderr, "                    the current time is used as seed.\n");
   fprintf(stderr, "  -q               Quiet mode\n");
   fprintf(stderr, "  -h               This message\n");

   if (ex>=0) exit(ex);
}

/* ============================================================ *
 * Main program.						*
 * ============================================================ */
int main(int argc, char *argv[])
{
   error *myerr = NULL, **err;
   config_pmc config;
   char *cname, pname[NSTR], iterdirname[NSTR];
   FILE *FLOG=NULL, *EVI=NULL, *PERP=NULL, *ENC=NULL, *PMCSIM, *PROP;
   time_t t_start;
   parabox *pb = NULL;
   gsl_rng *rng;
   int myid, nproc, iter, c, seed, this_nsamples=0, sum_nsamples=0, Nrevive=1, quiet;
   pmc_simu *psim;
   mix_mvdens *proposal;
   double sum_ess=0.0;
   extern char *optarg;
   extern int optind, optopt;


   /* ==================================================================== */

   err = &myerr;


   /* Command line options */
   cname = NULL;
   quiet = 0;
   seed  = -1;
   while ((c = getopt(argc, argv, ":c:s:qh")) != -1) {
      switch (c) {
         case 'c' :
            cname = optarg;
            break;
         case 's' :
            seed = atoi(optarg);
            break;
         case 'h' :
            usage(0, NULL);
         case 'q' :
            quiet = 1;
            break;
         case ':' :
            usage(1, "Argument for -%c option missing\n", optopt);
         case '?' :
            usage(2, "Unknown option -%c\n", optopt);
         default :
            /* Do nothing */
            break;
      }
   }
   if (cname==NULL) {
      cname = malloc_err(NSTR*sizeof(char), err);
      quitOnError(*err, __LINE__, stderr);
      strcpy(cname, "config_pmc");
   }

   /* ==================================================================== */


   /* Initialise  MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
   gsl_set_error_handler_off();

   sprintf(pname, "%s_pmc", log_name);
   FLOG = fopen_err(pname, "w", err); quitOnError(*err, __LINE__, stderr);

   /* Read configuration file */
   if (myid==0) {
      fprintf(FLOG, "%s v%s compiled on %s %s\n", __FILE__, COSMO_PMC_VERSION, __DATE__, __TIME__);
      time(&t_start);
      fprintf(FLOG, "started at %s\n", ctime(&t_start));
      fprintf(FLOG, "Number of nodes nproc = %d\n", nproc);
      if (! quiet) {
         fprintf(stderr, "%s v%s started on %d node%s\n", __FILE__, COSMO_PMC_VERSION, nproc, nproc>1?"s":"");
      }
   }

   if (myid!=0) {
      FLOG = fopen_err(pname, "a", err); quitOnError(*err, __LINE__, stderr);
   }


   rng = init_random(seed, FLOG);


   /* Alternative for reading config file without IO conflict */
   sleep(myid*0.01);
   /* Read config file and initialize proposal */
   read_config_pmc_file(&config, cname, &proposal, rng, 0, err);
   if (myid==0) {
      quitOnErrorStr(*err, __LINE__, stderr,
            getErrorValue(*err)==mv_cholesky ? 
            "***\nSuggestions:\n"
            "- Use the diagonal of the Fisher matrix (change 'sinitial' key in 'config_fish' to 'Fisher_diag'\n"
            "  or run 'cosmo_pmc.pl with option '-d')\n"
            "- Tweak the Fisher matrix by hand (e.g. decrease off-diagonal elements)\n***" : "");
   }

   if (myid==0) {
      fprintf(FLOG, "config file: %s\n", cname);
      out_config_pmc(FLOG, config, err);                         quitOnError(*err, __LINE__, stderr);

      pb = parabox_from_config(config.base.npar, config.base.min, config.base.max, err); quitOnError(*err, __LINE__, stderr);
   }



   /* ==================================================================== */


   /* New: Use _mpi init function for server and client */
   psim = pmc_simu_init_mpi(config.nsamples, config.base.npar, config.base.n_ded, err);
   quitOnError(*err, __LINE__, stderr);
   //fprintf(stderr, "myid = %d, psim = %d %d\n", myid, psim->mpi_rank, psim->mpi_size);


   /* Initialize PMC simulation, send and receive initial proposal */
   if (myid==0) {        /* Server */

      /* Evidence (Laplace approximation) */
      sprintf(pname, "%s_%s", evidence_name, fisher_suf);
      evidence_approx(&config.base, fisher_name, pname, err);
      checkPrintErr_and_continue(err, stderr, FLOG, NULL, 0);

      /* Open some files for writing */
      PERP = fopen_err(perplexity_name, "w", err); quitOnError(*err, __LINE__, stderr);
      EVI = fopen_err(evidence_name, "w", err);    quitOnError(*err, __LINE__, stderr);
      ENC = fopen_err(enc_name, "w", err);         quitOnError(*err, __LINE__, stderr);

      write_enc(proposal, 0, ENC, err);            quitOnError(*err, __LINE__, stderr);

      /* For perplexity/ess output */
      if (myid==0) sum_nsamples = sum_ess = 0;

   }


   /* ==================================================================== */


   /* Run PMC iterations */
   for (iter=0; iter<config.niter; iter++) {

      /* Directory for output files during the iteration */
      sprintf(iterdirname, "%s_%d", iter_dir, iter);
      make_and_test_dir(iterdirname, err);
      quitOnError(*err, __LINE__, stderr);

      if (myid==0) {
         this_nsamples = (int)(config.nsamples*(iter==config.niter-1 ? config.fsfinal : 1.0));
         sum_nsamples += this_nsamples;
      }

      /* Check wether files iter/pmcsim and iter/proposal exist */
      sprintf(cname, "%s/%s", iterdirname, pmcsim_name);
      sprintf(pname, "%s/%s", iterdirname, proposal_name);
      PMCSIM = fopen(cname, "r"); PROP   = fopen(pname, "r");

      if (PMCSIM==NULL || PROP==NULL) {  /* At least one file does not exist -> Run PMC iteration */

         if (PMCSIM) fclose(PMCSIM); if (PROP) fclose(PROP);

         run_pmc_iteration_MPI(psim, &proposal, iter, this_nsamples, &Nrevive, &config, iterdirname, myid, nproc,
               quiet, rng, pb, FLOG, err);
         quitOnError(*err, __LINE__, stderr);

      } else {                          /* Both PMC simulation and proposal files exist -> Read PMC iteration */

         if (psim) pmc_simu_free(&psim);
         psim = read_pmc_iteration_MPI(&proposal, iter, this_nsamples, &config, myid, nproc, PMCSIM, PROP, FLOG, err);
         quitOnError(*err, __LINE__, stderr);

         fclose(PMCSIM); fclose(PROP);

      }

      if (myid==0) {  /* Server */

         post_processing(psim, proposal, iter, sum_nsamples, &sum_ess, &config.base, iterdirname,
               PERP, EVI, ENC, FLOG, err);
         quitOnError(*err, __LINE__, stderr);

      }

   } /* iter */


   /* ==================================================================== */


   /* Wrap up */
   if (myid==0) {        /* Server */

      fclose(PERP);
      fclose(ENC);
      fclose(EVI);

      /* Final proposal */
      write_proposal(proposal, NULL, FLOG, -1, NULL, err); quitOnError(*err, __LINE__, stderr);
      free_parabox(&pb);

   }


   if (psim) pmc_simu_free(&psim);
   gsl_rng_free(rng);
   mix_mvdens_free(&proposal);     quitOnError(*err, __LINE__, stderr);


   if (myid==0) {       /* Server */
      fprintf(FLOG, "cosmo_pmc ");
      end_time(t_start, FLOG);
      fprintf(stderr, "cosmo_pmc finished\n");
   }

   fclose(FLOG);


   MPI_Finalize();

   return 0;
}


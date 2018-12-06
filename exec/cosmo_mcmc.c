/* ============================================================ *
 * cosmo_mcmc.c							*
 * Martin Kilbinger, Karim Benabed 2007/2008			*
 * ============================================================ */

#include "cosmo_mcmc.h"

/* ============================================================== *
 * One MCMC step. Returns 1 if new parameter vector is accepted.  *
 * On input logL is previous, on output it is new log-likelihood. *
 * ============================================================== */

#define NSTR 128
#define BASE 200

int mcmc_one(double *pstate, double *nstate, double *param_ded, double *logL, double *proba,
	     mc_law *mh, const parabox *pb, gsl_rng* rng,
	     config_base *config, FILE *FLOG, sm2_error **err)
{
  int j, inside;
  double transition, plogL, nlogL, r;
  MPI_Status status;
  
  plogL = *logL;
  
  /* Draw new parameter vector and transmits to slaves */
  memcpy(nstate, pstate, config->npar*sizeof(double));

  if(config->myid == 0){
    mh->sample(nstate, mh, rng, err);
    //print_parameter(stderr, config->npar, nstate);
    forwardError(*err, __LINE__, MC_AF_ALL);
    for(j=1;j<config->nproc;j++){
      MPI_Send(nstate, config->npar, MPI_DOUBLE, j, BASE+0, MPI_COMM_WORLD);
    }
  }else{
    MPI_Recv(nstate, config->npar, MPI_DOUBLE, 0, BASE+0, MPI_COMM_WORLD, &status);
  }
  

  inside = isinBox(pb, nstate, err);                          forwardError(*err, __LINE__, MC_AF_ALL);
  if (inside==0) {
    *logL = -MC_LOGBIG;
    return MC_AF_OUTOFBOX;
  }
  
  /* transition law */
  if (mh->transition_ratio!=NULL) {
    transition = mh->transition_ratio(pstate, nstate, mh, err);
    forwardError(*err, __LINE__, MC_AF_REJECT);
  } else {
    transition = 1;
  }
  
  /* Log-likelihood. NEW: can be computed in parallel */
  nlogL = posterior_log_pdf_common_MPI(config, nstate, err);
  //nlogL = posterior_log_pdf_common(config, nstate, err);
  
  *logL = nlogL;
  if (getErrorValue(*err)!=noErr) {
    fprintf(FLOG, "*** Error:\n");
    printError(FLOG, *err);
    purgeError(err);
    return MC_AF_REJECT;
  }
  
  forwardError(*err, __LINE__, MC_AF_ALL);
  
  /* transition probability */
  *proba = exp(nlogL-plogL)*transition;
  
  //fprintf(stderr, "proba (nlogL, plogL) = %g (%g %g)\n", *proba, nlogL, plogL);
  
   if (*proba>=1.0) {
     if (param_ded!=NULL && config->n_ded>0) {
       retrieve_ded((void*)config, param_ded, err);
       forwardError(*err, __LINE__, MC_AF_REJECT);
      }
     return MC_AF_ACCEPT;
   }
   
   /* the result is transmitted to guarantee identical parallel progression */
   /* (but it's maybe not necessary ? as only the file writing -- controlled by master -- is affected ) */
   if(config->myid == 0){
     r = gsl_ran_flat(rng, 0.0, 1.0);
     for(j=1;j<config->nproc;j++){
       MPI_Send(&r, 1, MPI_DOUBLE, j, BASE+1, MPI_COMM_WORLD);
     }
   }else{
     MPI_Recv(&r, 1, MPI_DOUBLE, 0, BASE+1, MPI_COMM_WORLD, &status);
   }   
   
   if (r<*proba) return MC_AF_ACCEPT;
   else return MC_AF_REJECT;
}
#undef BASE

#define BASE 300
double posterior_log_pdf_common_MPI(config_base *config, const double *x, error **err)
{
  double logl=0.0, logpost, logpost_id, logpr, pr, sigma_8=-1.0;
   double *xtmp;
   const double *xprior;
   int i, j;
   special_t special;
   common_like *like;
   MPI_Status status;

   for (i=0,logpost=0.0; i<config->ndata; i++) {
     /* one likehood per cpu or computed by master if ncpus < ndata */
     if(config->myid == i || (config->myid == 0 && i >= config->nproc)){
       
       /* Before lensing likelihood: get sigma_8 (from CAMB)  *
	* CMB has to come before lensing in the data set list */
       like = (common_like*)(config->data_extra[i]);
       if (config->data[i]==Lensing && sigma_8!=-1) ((lens_state*)(like->state))->sigma_8 = sigma_8;
       
	logl = likelihood_log_pdf_single(config->data_extra[i], x, config->data[i], err);
	forwardError(*err, __LINE__, 0);
	
	
	/* After CMB likelihood: save sigma_8 */
	if (config->data[i]==CMB) sigma_8 = ((wmap_state*)(like->state))->sigma_8;
	
	/* Prior. To be added for each data set. */
	if (config->data[i]!=Mring) {
	  logpr = config->logpr_default;
	} else {
	  logpr = 0.0;
	}
	
	/* Special priors */
	special = get_special(config->data[i], config->data_extra[i], err);
	forwardError(*err, __LINE__, 0.0);
	logpr += prior_log_pdf_special(special, config->par, config->min, config->max, config->npar, err);
	forwardError(*err, __LINE__, 0.0);
	
	logpost += logl + logpr;
     }
   }
   
   if(config->myid != 0){
     MPI_Send(&logpost, 1, MPI_DOUBLE, 0, BASE+0, MPI_COMM_WORLD);
   }else{
     for(j=1;j<config->nproc;j++){
       MPI_Recv(&logpost_id, 1, MPI_DOUBLE, j, BASE+0, MPI_COMM_WORLD, &status);
       logpost += logpost_id;
     }
   }
   
   /* Additional prior (e.g. previous experiment) */
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
   
   /* the result is transmitted to guarantee identical parallel progression */
   if(config->myid == 0){
     for(j=1;j<config->nproc;j++){
       MPI_Send(&logpost, 1, MPI_DOUBLE, j, BASE+1, MPI_COMM_WORLD);
     }
   }else{
     MPI_Recv(&logpost, 1, MPI_DOUBLE, 0, BASE+1, MPI_COMM_WORLD, &status);
   }
   
   MPI_Barrier(MPI_COMM_WORLD);
   return logpost;

}
#undef BASE


void print_matrix(const double *a, int n, FILE *F)
{
   int i, j;

   for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
	 fprintf(F, "% .3e ", a[i*n+j]);
      }
      fprintf(F, "\n");
   }
}

mvdens *initial_proposal(config_mcmc config, FILE *FLOG, FILE *CHPRE, double correct, error **err)
{
   mvdens *mhgauss, *fish;
   int i, j;
   FILE *FISH;
   double *accstates;
   long naccepted;
   double *pmean, *pvar;
   fpos_t headpos;

   if(config.base.myid == 0){
     fprintf(stderr, "Initial proposal: ");
   }
   
   mhgauss = mvdens_alloc(config.base.npar, err);    

   switch (config.base.initial) {

      case mcmcini_previous :

	 /* Previous chain */
	if(config.base.myid == 0){
	  testErrorRetVA(CHPRE==NULL, mcmc_file,
			 "Previous chain ('%s.%s') does not exist. See 'sinitial' key in the config file.",
			 *err, __LINE__, NULL, chainout, chainpre_suf);
	  fprintf(stderr, "Reading previous chain..."); fflush(stderr);
	}
	accstates = read_mkmc_chain(config.base.npar, 0, config.nchain, CHPRE, &naccepted, &headpos,
				    NULL, NULL, NULL, err);
	forwardError(*err, __LINE__, NULL);
	
	/* Rewind chain to after header */
	 fsetpos(CHPRE, &headpos);

	 // TODO: If nprevious<ncov, do we still use previous chain to get covariance?

	 pmean     = calloc_err(config.base.npar, sizeof(double), err);
	 forwardError(*err, __LINE__, NULL);
	 pvar      = calloc_err(config.base.npar*config.base.npar, sizeof(double), err);
	 forwardError(*err, __LINE__, NULL);

	 estimate_param_covar_weight(config.base.npar, naccepted, 0, accstates, NULL, pmean, pvar, err);
	 forwardError(*err, __LINE__, NULL);

	 free(accstates);
	 mvdens_from_meanvar(mhgauss, NULL, pvar, 1.0);
	 break;

      case mcmcini_diag :
	
	 /* Diagonal MV Gaussian with fixed variance */
	if(config.base.myid == 0){
	  fprintf(stderr, "Diagonal mvdens\n");
	}
	for (i=0; i<config.base.npar; i++) {
	  mhgauss->std[i*config.base.npar+i] = (config.base.max[i]-config.base.min[i])/config.boxdiv;
	}
	break;
	
      case mcmcini_fisher_diag :

	if(config.base.myid == 0){
	  fprintf(stderr, "Warning: sinitial = 'Fisher_diag' is treated as 'Fisher'\n");
	}
	  /* No break here! */

      case mcmcini_fisher : case mcmcini_fisher_inv :

	 FISH = fopen_err(fisher_name, "r", err);
	 forwardError(*err, __LINE__, NULL);
	 if(config.base.myid == 0){
	   fprintf(stderr, "Reading Fisher matrix\n");
	 }
	 fish = mvdens_dwnp(FISH, err);                 forwardError(*err, __LINE__, NULL);
	 if(config.base.myid == 0){
	   testErrorRetVA(fish->ndim!=config.base.npar, mk_npar,
			  "Fisher matrix has dimension %d, different from npar=%d (config file",
			  *err, __LINE__, NULL, fish->ndim, config.base.npar);
	   fprintf(FLOG, "Fisher matrix (from file)\n");
	   mvdens_print(FLOG, fish); fprintf(FLOG, "\n"); fflush(FLOG);
	 }
	 fclose(FISH);
	 
	 if (config.base.initial==mcmcini_fisher) {
	    sm2_inverse(fish->std, config.base.npar, err);
	    forwardError(*err, __LINE__, NULL);
	 }

	 if(config.base.myid == 0){
	   //	   mvdens_print(FLOG, fish);
	   fprintf(FLOG, "sqrt[F^{-1}_aa] = ");
	   for (i=0; i<config.base.npar; i++) fprintf(FLOG, "%g ", sqrt(fish->std[i*config.base.npar+i]));
	   fprintf(FLOG, "\n\n");
	 }
	 mvdens_from_meanvar(mhgauss, NULL, fish->std, 1.0);

	 mvdens_free(&fish);
	 break;

      default :
	if(config.base.myid == 0){
	  *err = addErrorVA(mcmc_unknown, "Invalid sinitial %d(%s) in config file", *err, __LINE__,
			    config.base.initial, config.base.sinitial);
	}
	return NULL;
	
   }


   for (i=0; i<config.base.npar; i++) {
      for (j=0; j<config.base.npar; j++) {
	 mhgauss->std[i*config.base.npar+j] = mhgauss->std[i*config.base.npar+j]*correct;
      }
   }

   if(config.base.myid == 0){
     fprintf(FLOG, "Initial proposal (multiplied with fudge factor %g:\n", correct);
     mvdens_print(FLOG, mhgauss); fprintf(FLOG, "\n"); fflush(FLOG);
   }

   mvdens_cholesky_decomp(mhgauss, err);             forwardError(*err, __LINE__, NULL);

   for (i=0; i<config.base.npar; i++) {
     if(config.base.myid == 0){
       testErrorRet(!finite(mhgauss->mean[i]), mcmc_infnan,
		    "inf or nan encouterned in Fisher matrix file (mean)",
		    *err, __LINE__, NULL);
     }
       for (i=0; i<config.base.npar; i++) {
	 for (j=0; j<config.base.npar; j++) {
	   if(config.base.myid == 0){
	     testErrorRet(!finite(mhgauss->std[i*config.base.npar+j]), mcmc_infnan, 
			  "inf or nan encouterned in Fisher matrix file (covariance)", *err, __LINE__, NULL);
	      }
	 }
       }
   }
   
   return mhgauss;
}


#define BASE 400

void run_mcmc(config_mcmc config, const parabox *pb, FILE *FLOG, FILE *CHPRE, int seed, error **err)
{
   int i, inside, reading_chpre, eof, j, nok;
   int lacc, nacc;			   /* Counter of accepted par. vectors */
   int accept=0, naccepted;
   int ibox, iboxmax=200;
   double *pstate, *nstate, *lstate;	   /* Previous, new, last accepted parameter vectors */
   double *param_ded;			   /* Deduced parameter vector */
   double *accstates;			   /* List of ncov previous accepted parameter vectors */
   double plogL=0.0, nlogL, llogL=0.0;     /* Previous, new, last accepted chi^2 */
   double *pmean, *pvar;		   /* Mean and variance of chain */
   double correct, proba, daccept;
   mvdens *mhgauss;                        /* Gaussian proposal */
   mc_law *mh;				   /* Proposal */
   FILE *CHALL, *CHACC;		           /* Chain output */
   FILE *F, *MEANVAR;
   char fname[100], dname[100];
   gsl_rng *rng;
   MPI_Status status;
   
 
   if (config.base.myid == 0) {
     fprintf(FLOG, "run_mcmc:\n");
     fprintf(FLOG, "=========\n");
     fflush(FLOG);
   }

   /* No gsl error reporting */
   gsl_set_error_handler_off();


   /* Fudge factor */
   correct = dsqr(config.fudge)/config.base.npar;

   /* Sampling law (proposal) = mv-Gaussian */
   mhgauss = initial_proposal(config, FLOG, CHPRE, correct, err);            forwardError(*err, __LINE__,);
   mh      = mc_mvdens_init(mhgauss, err);                                    forwardError(*err, __LINE__,);
   if (config.base.myid == 0) {
     fprintf(FLOG, "Initial proposal (after Cholesky):\n");
     mvdens_print(FLOG, mhgauss); fprintf(FLOG, "\n");
     fflush(FLOG);
   }
   pmean   = calloc_err(config.base.npar, sizeof(double), err);                    forwardError(*err, __LINE__,);
   pvar    = calloc_err(config.base.npar*config.base.npar, sizeof(double), err);   forwardError(*err, __LINE__,);




   /* Quit after Fisher matrix */
   if (config.start==start_nul) exit(0);
   if (config.nchain==0) exit(0);
   
   if (config.base.myid == 0) {
     fprintf(stderr, "Creating new chain\n");
   }
   
   accstates = calloc_err(config.base.npar*config.nchain, sizeof(double), err);
   forwardError(*err, __LINE__,);

   pstate    = accstates;
   nstate    = accstates+config.base.npar;
   lstate    = calloc_err(config.base.npar, sizeof(double), err);
   forwardError(*err, __LINE__,);

   param_ded = init_param_ded(config.base.n_ded, 1, err);
   forwardError(*err, __LINE__,);

   if (config.base.myid == 0) {
     /* ===== Open files for writing ===== */
     
     /* Chain output files */
     sprintf(fname, "%s/%s.%s", config.dirname, chainout, chainall_suf);
     CHALL = fopen_err(fname, "w", err); forwardError(*err,__LINE__,);
     
     sprintf(fname, "%s/%s.%s", config.dirname, chainout, chainacc_suf);
     CHACC = fopen_err(fname, "w", err);  forwardError(*err,__LINE__,);
     
     /* Mean and variance, updated at each covariance update */
     sprintf(fname, "%s/mean_update", config.dirname);
     MEANVAR = fopen_err(fname, "w",err );  forwardError(*err,__LINE__,);
     
   /* ================================== */
   
     write_header_cosmo_pmc(CHALL, config.base.par, config.base.npar, config.base.n_ded);
     write_header_cosmo_pmc(CHACC, config.base.par, config.base.npar, config.base.n_ded);
   }


 


   /* Initial parameter vector */
   if (config.base.myid == 0) {
     rng = init_random(seed, FLOG);
     
     while (!CHPRE) {
       set_start_parameter(pstate, config.start, config.fid, config.base, rng, FLOG, err);
       inside = isinBox(pb, pstate, err);            forwardError(*err, __LINE__,);
       testErrorRet(!inside, mcmc_outOfBound,
		    "Initial parameter not in box (min, max in config file). Check fid in config file",
		    *err, __LINE__,);
       /* Initial chi^2 */
       /* ATTENTION: computed only by master posterior_log_pdf_common and posterior_log_pdf_common_MPI
	* must lead to same result */
       plogL = posterior_log_pdf_common(&config.base, pstate, err);
       
       if (checkPrintErr_and_continue(err, FLOG, stderr, pstate, config.base.npar)) {
	 if (config.start==start_ran) fprintf(stderr, "Trying other random start parameter\n");
	 else {
	   fprintf(stderr, "Change the start parameter and run cosmo_mcmc again\n");
	   exit(5);
	 }
       } else {
	 break;
       }
     }
   }   
   
   /* send initial point to everyone */
   if(config.base.myid == 0){
     for(j=1;j<config.base.nproc;j++){
       MPI_Send(&plogL, 1, MPI_DOUBLE, j, BASE+0, MPI_COMM_WORLD);
       MPI_Send(pstate, config.base.npar, MPI_DOUBLE, j, BASE+1, MPI_COMM_WORLD);
     }
   }else{
     MPI_Recv(&plogL, 1, MPI_DOUBLE, 0, BASE+0, MPI_COMM_WORLD, &status);
     MPI_Recv(pstate, config.base.npar, MPI_DOUBLE, 0, BASE+1, MPI_COMM_WORLD, &status);
   }
   
   
   if (param_ded!=NULL && config.base.n_ded>0) {
      retrieve_ded((void*)&(config.base), param_ded, err);
      forwardError(*err, __LINE__,);
   }


   /* First point does not go into accept list */
   if (config.base.myid == 0) {
     print_step_cosmo_pmc(CHALL, 0.0, plogL, config.base.npar, config.base.n_ded, pstate, param_ded);
   }
   //goto end;

   nacc = 0;
   lacc = 1;
   naccepted = 0;

   if (CHPRE!=NULL && !feof(CHPRE)) {
      /* Overwrite initial parameter vector by first entry of previous chain */
     eof = read_next_mkmc_step(CHPRE, config.base.npar, config.base.n_ded, pstate, param_ded, &plogL, &daccept, err);
      forwardError(*err, __LINE__,);
      accept = (int)daccept;
      if (config.base.myid == 0) {
	print_step_cosmo_pmc(CHALL, accept, plogL, config.base.npar, config.base.n_ded, pstate, param_ded);
	print_step_cosmo_pmc(CHACC, accept, plogL, config.base.npar, config.base.n_ded, pstate, param_ded);
      }
      reading_chpre = 1;
   } else {
     reading_chpre = 0;
   }
   
   for (i=1; i<config.nchain; i++) {
     ibox = 0;
     do {
       if (reading_chpre==1) {
	 eof = read_next_mkmc_step(CHPRE, config.base.npar, config.base.n_ded, nstate, param_ded, &nlogL, &daccept, err);
	 forwardError(*err, __LINE__,);
	 accept = (int)daccept;
	 if (eof!=0) {
	   reading_chpre = 0;
	   purgeError(err);
	 }
	 forwardError(*err, __LINE__,);
       }
       if (reading_chpre==0) {
	 nlogL = plogL;
	 accept = mcmc_one(pstate, nstate, param_ded, &nlogL, &proba, mh, pb, rng, &config.base, FLOG, err);
	 checkPrintErr_and_continue(err, FLOG, stderr, pstate, config.base.npar);
       }
       
       testErrorRetVA(++ibox>=iboxmax, mk_boxmax, "The last %d proposed parameter outside of box at chain step #%d",
		      *err, __LINE__,, ibox, i);
       
     } while (accept==MC_AF_OUTOFBOX);
     
     
     if (accept==MC_AF_REJECT && naccepted>0) {
       /* Write last accepted parameter again to increase its weight */
       if (config.base.myid == 0) {
	 print_step_cosmo_pmc(CHACC, (double)MC_AF_ACCEPT, llogL, config.base.npar, config.base.n_ded, lstate, param_ded);
       }
     }
     
     if (config.base.myid == 0) {
       print_step_cosmo_pmc(CHALL, (double)accept, nlogL, config.base.npar, config.base.n_ded, nstate, param_ded);
     }
     
     
     
     if (accept==MC_AF_ACCEPT) {
       /* Move one up in acceptance list */
       pstate = nstate;
       memcpy(lstate, nstate, config.base.npar*sizeof(double));
       llogL  = nlogL;
       
       if (config.base.myid == 0) {
	 print_step_cosmo_pmc(CHACC, (double)accept, nlogL, config.base.npar, config.base.n_ded, nstate, param_ded);
       }
       
       lacc   = (lacc+1)%config.ncov;
       nstate = accstates+config.base.npar*lacc;
       plogL  = nlogL;
       naccepted++;
       
       /* Progress information to stderr */
       if (0 && naccepted%10==0 && config.base.myid == 0) {
	 fprintf(stderr, "\r%7d %7d %7.2f     \r", i, naccepted, naccepted/(double)i);
	 fflush(stderr);
       }
       
       //if (naccepted==1) goto end;
       
       /* (re)calculate covariance */
       if (lacc==0) {

	 //fprintf(stderr, "naccepted = %d\n", naccepted);
	    estimate_param_covar_weight(config.base.npar, naccepted, config.ncov*nacc, accstates, NULL, pmean, pvar, err);
	    forwardError(*err, __LINE__,);
	    nacc++;
	    
	    if (config.base.myid == 0) {
	      fprintf(MEANVAR, "%5d", naccepted);
	      for (j=0; j<config.base.npar; j++) {
		fprintf(MEANVAR, "   % .5g % .5g", pmean[j], pvar[j*config.base.npar+j]);
	      }
	      fprintf(MEANVAR, "\n"); fflush(MEANVAR);
	      fprintf(stderr, "new cov: %f %f\n", pvar[0], correct);
	    }
	    
	    mvdens_from_meanvar(mhgauss, NULL, pvar, correct);
	    
	    if (config.base.myid == 0) {
	      fprintf(stderr, "mvdens: %f %f\n", mhgauss->std[0], correct);
	      fprintf(FLOG, "step %d, acc %d,  parameter covariance (multiplied with %g):\n", i, naccepted, correct);
	      mvdens_print(FLOG, mhgauss); fprintf(FLOG, "\n");
	      
	      if (OUT_COV==1) {
		sprintf(dname, "%s/%s", config.dirname, covar_dir);
		mkdir(dname, 0755);
		nok = chdir(dname);
		if (nok==0) {    /* chdir successful */
		  sprintf(fname, "cov_%06d", naccepted);
		  F = fopen_err(fname, "w", err);                    forwardError(*err, __LINE__,);
		  mvdens_dump(F, mhgauss);
		  fclose(F);
		  
		  sm2_inverse(mhgauss->std, config.base.npar, err);       forwardError(*err, __LINE__,);
		  sprintf(fname, "covinv_%06d", naccepted);
		  F = fopen_err(fname, "w", err);                    forwardError(*err, __LINE__,);
		  mvdens_dump(F, mhgauss);            
		  fclose(F);
		  sm2_inverse(mhgauss->std, config.base.npar, err);       forwardError(*err, __LINE__,);
		  
		  chdir("../..");
		} else {
		  fprintf(stderr, "Warning: Could not create and change into directory %s. Make sure you have "
			  "write permission.\n", fname);
		}
	      }
	      fprintf(stderr, "after out: %f %f\n", mhgauss->std[0], correct);
	    }
	    
	    mvdens_cholesky_decomp(mhgauss, err);
	    forwardError(*err, __LINE__,);
	    if (config.base.myid == 0) {
	      fprintf(stderr, "Step %d, acc %d,  parameter covariance (multiplied with %g):\n", i, naccepted, correct);
	      mvdens_print(stderr, mhgauss);
	      fprintf(FLOG, "Step %d, acc %d,  parameter covariance (after cholesky):\n", i, naccepted);
	      mvdens_print(FLOG, mhgauss); fprintf(FLOG, "\n");
	    }
       }
     }
   }
   
   // end:
   
   if (config.base.myid == 0) {
     fclose(CHALL);
     fclose(CHACC);
     fclose(MEANVAR);
     gsl_rng_free(rng);
   }
   
   
   if (config.base.myid == 0) {
     fprintf(FLOG, "naccepted = %d, rate = %.3f\n", naccepted, (double)(naccepted)/config.nchain); fflush(FLOG);
     fprintf(stderr, "naccepted = %d, rate = %.3f\n", naccepted, (double)(naccepted)/config.nchain); fflush(stderr);
     sprintf(fname, "%s/rate", config.dirname);
     F = fopen_err(fname, "w",err);  forwardError(*err,__LINE__,);
     fprintf(F, "naccepted = %d, rate = %.3f\n", naccepted, (double)(naccepted)/config.nchain);
     fclose(F);
     fprintf(FLOG, "\n");
   }
   
   free(pmean);
   free(pvar);
   free(lstate);
   mc_mvdens_free(&mh);
   free(accstates);
}

#undef BASE

double *burn_in_decorrelate(double *accstates, double *accparam_ded, long naccepted, long *nfinal,
			    config_mcmc config, FILE *FLOG, double **finalparam_ded, error **err)
{
   FILE *CHFIN;
   double *finalstates;
   long iacc, ifin;
   char fname[512];

   testErrorRet(config.ndecorr==0, mcmc_negative,
		"ndecorr (in config_mcmc) has to be one or larger", *err, __LINE__, NULL);

   fprintf(stderr, "Removing burn-in and decorrelating\n");
   fprintf(FLOG, "burn_in_decorrelate:\n");
   fprintf(FLOG, "====================\n");

   finalstates    = calloc_err(config.base.npar*naccepted, sizeof(double), err);
   forwardError(*err, __LINE__, 0);
   if (accparam_ded!=NULL) {
      *finalparam_ded = calloc_err(config.base.n_ded*naccepted, sizeof(double), err);
      forwardError(*err, __LINE__, 0);
   } else {
      *finalparam_ded = NULL;
   }

   sprintf(fname, "%s/%s.%s", config.dirname, chainout, chainfin_suf);
   CHFIN = fopen_err(fname, "w",err);  forwardError(*err,__LINE__,NULL);
   write_header_cosmo_pmc(CHFIN, config.base.par, config.base.npar, config.base.n_ded);

   ifin = 0;
   iacc = (int)(config.ncov*config.fburnin);
   fprintf(FLOG, "Discarding first %ld points\n", iacc);

   while (iacc<naccepted) {
      memcpy(finalstates+ifin*config.base.npar, accstates+iacc*config.base.npar, config.base.npar*sizeof(double));
      if (accparam_ded!=NULL) {
	 memcpy((*finalparam_ded)+ifin*config.base.n_ded, accparam_ded+iacc*config.base.n_ded, config.base.n_ded*sizeof(double));
      }

      print_step_cosmo_pmc(CHFIN, (double)MC_AF_ACCEPT, 0.0, config.base.npar, config.base.n_ded, 
			   accstates+iacc*config.base.npar, accparam_ded+iacc*config.base.n_ded);
      iacc += config.ndecorr;
      ifin ++;
   }

   fprintf(stderr, "After decorrelation %ld points left\n\n", iacc);

   fclose(CHFIN);
   fprintf(FLOG, "After decorrelation %ld points left\n\n", iacc);
   *nfinal = ifin;

   return finalstates;
}

/* Mean and covariance as fct of time */
void mean_cov_time(double *states, int nstates, int every, config_base config, error **err)
{
   FILE *F, *F2;
   int n, i, j;
   double *pmean, *pvar;

   pmean   = calloc_err(config.npar, sizeof(double), err);                forwardError(*err, __LINE__,);
   pvar    = calloc_err(config.npar*config.npar, sizeof(double), err);    forwardError(*err, __LINE__,);
   F  = fopen_err("meanT", "w",err);  forwardError(*err,__LINE__,);
   F2 = fopen_err("varT", "w",err);  forwardError(*err,__LINE__,);

   for (n=every; n<nstates; n+=every) {
      for (i=0; i<config.npar; i++) {
	 pmean[i] = 0.0;
	 for (j=0; j<config.npar; j++) {      
	    pvar[i*config.npar+j] = 0.0;
	 }
      }
      estimate_param_covar_weight(config.npar, n, 0, states, NULL, pmean, pvar, err);
      forwardError(*err, __LINE__,);

      fprintf(F, "%5d", n);
      for (i=0; i<config.npar; i++) {
	 fprintf(F, " %.7g", pmean[i]);
      }
      fprintf(F, "\n");

      fprintf(F2, "%5d", n);
      for (i=0; i<config.npar; i++) {
	 fprintf(F2, " %.7g", pvar[i*config.npar+i]);
      }
      fprintf(F2, "\n");
   }

   fclose(F);
   fclose(F2);
   free(pmean);
   free(pvar);
}

double mean_from_chain(const double *states, long nstates, int npar, int a)
{
   int i;
   double m;

   for (i=0,m=0.0; i<nstates; i++) {
      m += states[i*npar + a];
   }
   m /= (double)nstates;

   return m;
}

void sigma_from_chain(const double *states, int nstates, int npar, int a, double mean, double *sigma, error **err)
{
   const double confidence[3] = {0.6827, 0.9545, 0.9973}; 
   int i, j, imean;
   double sum, invn;
   double *tmp;

   /* Copy relevant data and sort according to parameter values */
   tmp = malloc_err(sizeof(double)*nstates, err);   forwardError(*err, __LINE__,);
   for (i=0; i<nstates; i++) {
      tmp[i] = states[i*npar+a];
   }
   qsort(tmp, nstates, sizeof(double), double_cmp);

   /* Look for mean */
   i = 0;
   while (tmp[i]<mean) {
      testErrorRet(i==nstates-1, mk_range, "Mean not found, maybe all chain points are the same?",
		   *err, __LINE__,);
      i++;
   }
   imean = i;

   /* Sum up right of mean until confidence volume */
   invn = 1.0/(double)nstates;
   for (j=0; j<3; j++) {
      sum = 0.0;
      i = imean;
      while (sum<=confidence[j]/2.0 && i<nstates-1) {
	 sum += invn;
	 //fprintf(stderr, "%d %g %g %g  %g\n", i, invn, sum, confidence[j]/2.0, mean-tmp[i]);
         i++;
      }
      sigma[j] = tmp[i]-mean;
      if (i==nstates-1) sigma[j] = -sigma[j];
      //fprintf(stderr, "sigma = %g\n", sigma[j+3]);
   }

   /* Sum up left of mean until confidence volume */
   for (j=0; j<3; j++) {
      sum = 0.0;
      i = imean-1;
      while (sum<=confidence[j]/2.0 && i>0) {
	 sum += invn;
         i--;
      }
      sigma[j+3] = mean-tmp[i];
      if (i==0) sigma[j+3] = -sigma[j+3];
   }

   free(tmp);
}

void mean_sigma_from_chain(config_mcmc config, const double *states, const double *states_ded, long nstates, int npar, int n_ded,
			   error **err)
{
   FILE *F;
   double *pmean, *pvar;
   int i, j;
   char fname[NSTR];

   sprintf(fname, "%s/%s", config.dirname, mean_name);
   
   F = fopen_err(fname, "w",err);  forwardError(*err,__LINE__,);
   pmean = calloc_err(npar, sizeof(double), err);          forwardError(*err, __LINE__,);
   pvar  = calloc_err(npar*6, sizeof(double), err);        forwardError(*err, __LINE__,);

   for (i=0; i<npar; i++) {
      pmean[i] = mean_from_chain(states, nstates, npar, i);
      sigma_from_chain(states, nstates, npar, i, pmean[i], pvar, err);
      forwardError(*err, __LINE__,);
      fprintf(F, "%2d %10.5g", i, pmean[i]);
      for (j=0; j<3; j++) fprintf(F, "   % 9.5g % 9.5g", pvar[j], pvar[j+3]);
      fprintf(F, "\n");
   }

   for (i=0; i<n_ded; i++) {
      pmean[i] = mean_from_chain(states_ded, nstates, n_ded, i);
      sigma_from_chain(states_ded, nstates, n_ded, i, pmean[i], pvar, err);
      forwardError(*err, __LINE__,);
      fprintf(F, "%2d %10.5g", i+npar, pmean[i]);
      for (j=0; j<3; j++) fprintf(F, "   % 9.5g % 9.5g", pvar[j], pvar[j+3]);
      fprintf(F, "\n");
      }

   fclose(F);
   free(pmean);
   free(pvar);
}

double *init_param_ded(int n_ded, int nmaxchain, error **err)
{
   double *param_ded;

   if (n_ded>0) {
      param_ded = malloc_err(n_ded*nmaxchain*sizeof(double), err);
      forwardError(*err, __LINE__, NULL);
   } else {
      param_ded = NULL;
   }

   return param_ded;
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

   fprintf(stderr, "Usage: cosmo_mcmc [OPTIONS]\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: config_mcmc)\n");
   fprintf(stderr, "  -s SEED          Use SEED for random number generator. If SEED=-1 (default)\n");
   fprintf(stderr, "                    the current time is used as seed.\n");
   fprintf(stderr, "  -o DIR           Output directory (default: current)\n");
   fprintf(stderr, "  -h               This message\n");

   if (ex>=0) exit(ex);
}

/* ============================================================ *
 * Main program.						*
 * ============================================================ */


#define BASE 500
int main(int argc, char *argv[])
{
   error *myerr = NULL, **err;
   config_mcmc config;
   FILE *FLOG, *CHAIN, *CHPRE;
   char fname[NSTR], *cname;
   time_t t_start;
   double *accstates, *finalstates, *param_ded, *finalparam_ded;
   double *finalall;
   parabox *pb;
   long naccepted, nfinal;
   fpos_t headpos;
   int c, seed;
   extern char *optarg;
   extern int optind, optopt;

  /* Initialise  MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &config.base.myid);
   MPI_Comm_size(MPI_COMM_WORLD, &config.base.nproc);
   gsl_set_error_handler_off();

   /* ==================================================================== */
   
   err = &myerr;
   
   /* Command line options */
   cname = NULL;
   config.dirname = NULL;
   seed  = -1;
   while ((c = getopt(argc, argv, ":c:o:s:h")) != -1) {

      switch (c) {
      case 'c' :
	cname = optarg;
	break;
      case 'o' :
   	config.dirname = optarg;
	break;
      case 's' :
	seed = atoi(optarg);
	break;
      case 'h' :
	usage(0, NULL);      
      case ':' :
	usage(1, "Argument for -%c option missing\n", optopt);
      case '?' :
	usage(2, "Unknown option -%c\n\n", optopt);
      default :
	/* Do nothing */
	    break;
      }
   }

   if (cname==NULL) {
      cname = malloc_err(NSTR*sizeof(char), err);
      quitOnError(*err, __LINE__, stderr);
      strcpy(cname, "config_mcmc");
   }

   if (config.dirname==NULL) {
      config.dirname = malloc_err(NSTR*sizeof(char), err);
      quitOnError(*err, __LINE__, stderr);
      strcpy(config.dirname, ".");
   }
   
   /* ==================================================================== */

   if (config.base.myid==0) {
     fprintf(stderr, "cosmo_mcmc started\n");
   }



   
   /* Directory for output files */
   make_and_test_dir(config.dirname, err);
   quitOnError(*err, __LINE__, stderr);

   sprintf(fname, "%s/%s_mcmc", config.dirname, log_name);

   /* Alternative for reading config file without IO conflict */
   sleep(config.base.myid*0.01);
   read_config_mcmc_file(&config, cname, err);      quitOnError(*err, __LINE__, stderr);
   
   /* Log file */
   if (config.base.myid==0) {
     FLOG = fopen_err(fname, "w", err);  quitOnError(*err, __LINE__, stderr);
     fprintf(FLOG, "%s v%s compiled on %s %s\n", __FILE__, COSMO_PMC_VERSION, __DATE__, __TIME__);
     t_start = start_time(FLOG);
     
     fprintf(stderr, "Config file: %s\n", cname);
     fprintf(FLOG, "Config file: %s\n", cname);
     out_config_mcmc(FLOG, config, err);              quitOnError(*err, __LINE__, stderr);
     fflush(FLOG);
   }
   if (config.base.myid!=0) {
     FLOG = fopen_err(fname, "a", err); quitOnError(*err, __LINE__, stderr);
   }

   pb = parabox_from_config(config.base.npar, config.base.min, config.base.max, err); quitOnError(*err, __LINE__, stderr);

   sprintf(fname, "%s/%s.%s", config.dirname, chainout, chainfin_suf);
   CHAIN = fopen(fname, "r");
   
   
   if (!CHAIN) {
     
     /* chain.fin does not exist: check for chain.acc */
     sprintf(fname, "%s/%s.%s", config.dirname, chainout, chainacc_suf);
     
     CHAIN = fopen(fname, "r");
     if (!CHAIN) {
       /* chain.acc does not exist: run mcmc, create chain.acc */
       sprintf(fname, "%s/%s.%s", config.dirname, chainout, chainpre_suf);

       
 
	
       CHPRE = fopen(fname, "r");
       if (config.base.myid==0) {
	 if (CHPRE!=NULL) {
	   fprintf(stderr, "Continuing with previous chain\n");
	   fprintf(FLOG, "Continuing with previous chain\n");
	 }
       }





       /* wait for everyone before moving on (prevents stoping the process if 
	* master creates chain.acc and a slave skips this condition */
       MPI_Barrier(MPI_COMM_WORLD);
       
       /* main loop */
       run_mcmc(config, pb, FLOG, CHPRE, seed, err);  quitOnError(*err, __LINE__, stderr);
       
       if (CHPRE) fclose(CHPRE);
       
     } else {
       fclose(CHAIN);
     }
     
     if (config.base.myid==0) {

       /* Read chain.acc */
       sprintf(fname, "%s/%s.%s", config.dirname, chainout, chainacc_suf);
       fprintf(stderr, "Reading %s\n", fname);
       CHAIN = fopen_err(fname, "r",err);  quitOnError(*err, __LINE__, stderr);
       param_ded = init_param_ded(config.base.n_ded, config.nchain, err);
       accstates = read_mkmc_chain(config.base.npar, config.base.n_ded, config.nchain, CHAIN, &naccepted, &headpos,
				   param_ded, NULL, NULL, err);
       quitOnError(*err, __LINE__, stderr);
       fclose(CHAIN);
       
       testErrorRetVA(naccepted==0, mk_empty, "Number of accepted points is zero. Maybe chain '%s' is empty or only contains header?",
		      *err, __LINE__, 3, fname);
       
       fprintf(stderr, "Mean and cov as fct of time..."); fflush(stderr);
       //mean_cov_time(accstates, naccepted, config.nchain/1000, config, err);
       fprintf(stderr, "done\n");
       quitOnError(*err, __LINE__, stderr);
       
       /* Remove burn-in phase, decorrelate (thin-out) chain and create chain.fin */
       finalstates = burn_in_decorrelate(accstates, param_ded, naccepted, &nfinal, config, FLOG, &finalparam_ded, err);
       quitOnError(*err, __LINE__, stderr);
       free(accstates);
     }
     
     
   } else {
     
     if (config.base.myid==0) {
       
       /* chain.fin exists */
       fprintf(stderr, "Reading %s/%s.%s\n", config.dirname, chainout, chainfin_suf);
       param_ded = init_param_ded(config.base.n_ded, config.nchain, err);
       finalstates = read_mkmc_chain(config.base.npar, config.base.n_ded, config.nchain, CHAIN, &nfinal, &headpos,
				     param_ded, NULL, NULL, err);
       quitOnError(*err, __LINE__, stderr);
       /* Is not needed anymore, buy anyway, deal with deduced parameters correctly */
       finalparam_ded = param_ded;
       fclose(CHAIN);
     }
   }
   
   if (config.base.myid == 0) {
     
     /* Parameter means and variances */
     mean_sigma_from_chain(config, finalstates, finalparam_ded, nfinal, config.base.npar, config.base.n_ded, err);
     quitOnError(*err, __LINE__, stderr);
     
     if (nfinal==0) {
       *err = addErrorVA(mk_empty, "Final chain (%s.%s) is empty", *err, __LINE__, chainout, chainfin_suf);
       quitOnError(*err, __LINE__, stderr);
     }
     
     /* Histograms */
     fprintf(stderr, "nfinal = %ld   naccepted = %ld\n", nfinal, naccepted);
     sprintf(fname, "rm -f %s_*", hist_name);
     system(fname);
     
     finalall = merge_param_param_ded(finalstates, finalparam_ded, nfinal, config.base.npar, config.base.n_ded, err);
     quitOnError(*err, __LINE__, stderr);
     
     create_1d_2d_hist(finalall, nfinal, config.base.nbinhist, config.base.npar+config.base.n_ded, config.base.min,
		       config.base.max, NULL, FLOG, config.dirname, err); 
     /* New: no exitOnError */

     /* Covariance and inverse of final chain */
     covariance_from_sample_all(finalstates, NULL, config.base.npar, nfinal, 
				finalall, config.base.npar+config.base.n_ded,  config.dirname, err);
     quitOnError(*err, __LINE__, stderr);
     
     
     
     /* ==================================================================== */
     
     free(finalstates);
     if (config.base.n_ded>0) free(finalall);
     free(param_ded);
     if (finalparam_ded==NULL) free(finalparam_ded);
     free_parabox(&pb);
     
     //free_config(config, *err);
     //quitOnError(*err, __LINE__, stderr);
     
     fprintf(FLOG, "cosmo_mcmc ");
     end_time(t_start, FLOG);
     fclose(FLOG);
     fprintf(stderr, "cosmo_mcmc finished\n");
     
   }



   

   
   MPI_Finalize();
   
   
   return 0;
}
#undef BASE

/*
 *  other_mcmc_tst.c
 *  montecarlo
 *
 *  Created by Karim Benabed on 12/10/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "montecarlo.h"
#include "parabox.h"
#include "target_circ.h"

#include <stdio.h>

#define T_NDIM 10

int  main(int argc, char** argv) {
	
  gsl_rng *r;
  mix_mvrnorm *model;
  mvrnorm *mhgauss;
  mc_law *mh;
  long i,l,j,indim;
  double *pstate,*nstate,*exc;
  double pvalue,nvalue;
  parabox* pb;
  long accept,naccept,lacc,nacc;
  double *pmean,*pvar,*acc_params;
  double correct;
	
  int ndim=T_NDIM;
  int ntry=100000;
  int ntry1=2000;
	
  correct=2.4*2.4/ndim;
	
  error *err;
  err=NULL;
	
  FILE* resFile;
	
  // nuking gsl error reporting
  gsl_set_error_handler_off();
	
  /* Initialize the random seed (beware: different calls will result	*/
  /* in idnetical runs)							*/
  r = gsl_rng_alloc(gsl_rng_default);
	
  /* my model ! */
  model = init_posterior(ndim,&err);
  exitOnError(err,stderr);
  fprintf(stderr,"Model\n");
  mix_mvrnorm_print(stderr,model);
	
	
  /* MH law */
  pmean=calloc(ndim,sizeof(double));
  pvar=calloc(ndim*ndim,sizeof(double));
	
  mhgauss=mvrnorm_alloc(ndim,&err);
  exitOnError(err,stderr);
	
  for (i = 0; i < ndim; i++) {
    mhgauss->std[i*ndim+i]=1.;
  }
  mvrnorm_print(stderr,mhgauss);
  mvrnorm_cholesky_decomp(mhgauss,&err);
  exitOnError(err,stderr);
  fprintf(stderr,"Test\n");
	
  mh=mc_mvrnorm_init(mhgauss,&err);
  exitOnError(err,stderr);
	
  /* a box */
  pb = init_parabox(ndim,&err);
  exitOnError(err,stderr);
	
  for (i=0;i<ndim;i++) {
    add_slab(pb, i, -10, 20 ,&err); 
    exitOnError(err,stderr);
  }
  
	
  /* some inits */
  acc_params = calloc(ndim*ntry1,sizeof(double));
  lacc=0;
  nacc=0;
  pstate = acc_params;
  nstate = acc_params+ndim;
	
  /* init point */
  while (1) {
    // replace that by a random init point
    for (i=0;i<ndim;i++) {
      pstate[i]=3;
    }
    if (isinBox(pb,pstate,&err)) {
      break;
    }
    exitOnError(err,stderr);
  }
  fprintf(stderr,"ici\n");
	
  // a file for the accepted points
  resFile=fopen("resRunOther.dat","w");
  
	
  pvalue = posterior_log_pdf(model,pstate, &err);
  exitOnError(err,stderr);
		
  naccept=0;
  l=0;
	
  print_step(resFile,1,pvalue,ndim,pstate);
  
	
  for (l=1;l<ntry;l++) {
  
    nvalue=pvalue;
    accept = mcmc_step(pstate,nstate,&nvalue,posterior_log_pdf,model,mh,pb,r,&err);
   
    //print_step(stderr,accept,nvalue,ndim,nstate);
    print_step(resFile,accept,nvalue,ndim,nstate);
    exitOnError(err,stderr);
    if ((accept!=0) ) {
      pstate=nstate;
      lacc=(lacc+1)%ntry1;
      nstate=acc_params+ndim*lacc;
      pvalue=nvalue;
      naccept++;
      if ((lacc==0)) {
	/* estimate the covariance matrix */
	estimate_param_covar(ndim,naccept,ntry1*nacc,acc_params,pmean,pvar,&err);
	exitOnError(err,stderr);
	nacc++;
	printf("step %d -> ",l);
	for(i=0;i<ndim;i++) {
	  indim=i*ndim;
	  printf("%g ",pmean[i]);
	  for (j=0;j<ndim;j++) {
	    mhgauss->std[indim+j]=pvar[indim+j]*correct;
	  }
	}
	printf("\n");
	mvrnorm_print(stderr,mhgauss);
	mvrnorm_cholesky_decomp(mhgauss,&err);
      }
    }
		
  }
	
  printf("naccept = %ld\n",naccept);
  fclose(resFile);			
  free_posterior(&model);
  free_parabox(&pb);
  mc_mvrnorm_free(&mh);
	
}


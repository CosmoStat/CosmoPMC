/* ============================================================ *
 * xcorr_post.c							*
 * MK 2008							*
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
#include "parabox.h"
#include "param.h"
#include "config.h"
#include "mvdens.h"

#include "mcmc.h"
#include "cosmo_mcmc.h"
#include "stdnames.h"

#include "xcorr.h"


void transform_chain(config_mcmc config, const double *finalstates, int nfinal, error **err)
{
   int i, j, k, m;
   xcorr_state *xcorr;
   xcorr_model *model;
   double *params_trans, *N;
   FILE *CHTRA;
   char fname[100];

   xcorr = (xcorr_state*)config.data_extra[0];
   xcorr_print(stderr, xcorr);

   sprintf(fname, "%s.tra", chainout);
   CHTRA = fopen(fname, "w");
   testErrorRet(!CHTRA, mcmc_io, "Could not create chain output file", *err, __LINE__,);
   write_header_mk(CHTRA, xcorr->Mz*(xcorr->Mz+1));

   params_trans = vector(0, xcorr->Mz*(xcorr->Mz+1)-1, err);   forwardError(*err, __LINE__,);
   N            = vector(0, xcorr->Mz, err);	               forwardError(*err, __LINE__,);

   for (i=0; i<nfinal; i++) {
      model = params_to_xcorr_model(xcorr, finalstates+i*config.npar, err);
      forwardError(*err, __LINE__,);

      /* True redshift distribution N(z_j) */
      for (j=0; j<xcorr->Mz; j++) {
	 for (k=0,N[j]=0.0; k<xcorr->Mz; k++) {
	    N[j] += model->etainv[j][k]*xcorr->nbar[k];
	 }
      }

      for (j=0,m=0; j<xcorr->Mz; j++) {
	 for (k=0; k<xcorr->Mz; k++) {
	    //if (k==j) continue;
	    params_trans[m++] = model->etainv[k][j];
	 }
      }
      for (j=0; j<xcorr->Mz; j++) {
	 params_trans[m++] = N[j];
      }

      if (model) free_xcorr_model(model);
      print_step_mk(CHTRA, MC_AF_ACCEPT, 0.0, xcorr->Mz*(xcorr->Mz+1), params_trans);
   }


   free_vector(params_trans);
   free_vector(N);
   free(xcorr);

   fclose(CHTRA);
}

int main(int argc, char *argv[])
{
   error *err = NULL;
   config_mcmc config;
   FILE *CHAIN;
   char fname[100];
   double  *finalstates;
   int nfinal;
   fpos_t headpos;


   fprintf(stderr, "%s compiled on %s %s\n", __FILE__, __DATE__, __TIME__);


   read_config_file(&config, "config", &err); exitOnError(err, stderr);
   if (config.data[0]!=Xcorr) {
      fprintf(stderr, "data set has to be Xcorr");
      exit(2);
   }

   sprintf(fname, "%s.%s", chainout, chainfin_suf);
   CHAIN = fopen(fname, "r");
   if (!CHAIN) {
      fprintf(stderr, "%s.%s not found", chainout, chainfin_suf);
      exit(mcmc_io);
   }

   finalstates = read_mkmc_chain(config.npar, config.nchain, CHAIN, stderr, &nfinal, &headpos, &err);
   exitOnError(err, stderr);
   fclose(CHAIN);


   transform_chain(config, finalstates, nfinal, &err); exitOnError(err, stderr);


   return 0;
}


/* ============================================================ *
 * wmap5_post.c							*
 * Martin Kilbinger 2008					*
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
#include "config.h"
#include "mvdens.h"

#include "param.h"
#include "mcmc.h"
#include "cosmo_mcmc.h"
#include "stdnames.h"

#include "wmap.h"
#include "lensing.h"
#include "C_wrappers.h"

#define eps 1.0e-6
void transform_chain(config_mcmc config, const double *finalstates, int nfinal, error **err)
{
   int i, Ntrans, ic, j, same;
   wmap_like *wmap;
   forked_state *lkl;
   wmap_state *state;
   double *params_trans, *last;
   FILE *CHTRA;
   char fname[100];

   Ntrans = 6;

   fprintf(stderr, "Initialising wmap...\n");
   wmap = wmap_init("scamb", "data/", config.param, TT_MAX,0,NULL,DEFAULT_TIMEOUT, err);
   fprintf(stderr, "done\n");
   lkl = wmap->lkl;
   state = (wmap_state*)lkl->specific_data;
   wmap->want_sig8 = 1;
   wmap->want_out_Cl = 0;

   sprintf(fname, "%s.tra", chainout);
   CHTRA = fopen(fname, "w");
   testErrorRet(!CHTRA, mcmc_io, "Could not create chain output file", *err, __LINE__,);
   //write_header_lc(CHTRA, Ntrans);

   params_trans = calloc(Ntrans, sizeof(double));
   testErrorRet(params_trans==NULL, mcmc_allocate, "Could not allocate memory", *err, __LINE__,);
   last         = calloc(config.npar, sizeof(double));
   testErrorRet(last==NULL, mcmc_allocate, "Could not allocate memory", *err, __LINE__,);

   for (i=0; i<nfinal; i++) {

      ic = i*config.npar;

      same = 1;
      for (j=0; j<config.npar; j++) {
	 if (fabs(finalstates[ic+j]-last[j])>eps) {
	    same = 0;
            break;
      	 }
      }

      switch (config.param) {

	 case param_wCDM_flat_nphys :

	    /* 0  1  2   3  4  5  6 *
	     * Ob Oc tau w0 ns As h */

	    if (same==0) {

	       /* Omegab */
	       params_trans[0] = finalstates[ic+0];
	       /* Omegam = Omegab + Omegac */
	       params_trans[1] = finalstates[ic+0]+finalstates[ic+1];
	       params_trans[2] = finalstates[ic+2];
	       params_trans[3] = finalstates[ic+3];
	       params_trans[4] = finalstates[ic+4];
	       params_trans[6] = finalstates[ic+6];

	       /* sigma_8 */
	       fprintf(stderr, "Calling wmap_likely (want_sig8=1)...\n");
	       params_trans[5] = wmap_likely(wmap, finalstates+ic, err);
	       forwardError(*err, __LINE__,);
	       fprintf(stderr, "done\n");
	       fprintf(stderr, "sigma_8 = %g\n", params_trans[5]);
	    }
	    break;

         case param_WMAP5 :
	   /* 0     1  2 3  4   5
	    * 100wb oc h ns tau As */
	    
	   if (same==0) {
	      /* Omegab */
	      params_trans[0] = finalstates[ic+0]/100.0/dsqr(finalstates[ic+2]);
	      /* Omegam */
	      params_trans[1] = params_trans[0] + finalstates[ic+1]/dsqr(finalstates[ic+2]);
	      /* h */
	      params_trans[2] = finalstates[ic+2];
	      params_trans[3] = finalstates[ic+3];
	      params_trans[4] = finalstates[ic+4];
	      /* sigma_8 */
	      params_trans[5] = wmap_likely(wmap, finalstates+ic, err);
	      forwardError(*err, __LINE__,);
	      fprintf(stderr, "sigma_8 = %g\n", params_trans[5]);
	   }
	   break;

         case param_WMAP5_flat :
	   /* 0  1 2  3   4  5
	    * oc h ns tau As 100ob */
	    
	   if (same==0) {
	      /* Omegam = Omegab + omegac/h2 */
	      params_trans[0] = params_trans[0] + finalstates[ic+0]/dsqr(finalstates[ic+1]);
	      /* h */
	      params_trans[1] = finalstates[ic+1];
	      /* n_s */
	      params_trans[2] = finalstates[ic+2];
	      /* tau */
	      params_trans[3] = finalstates[ic+3];
	      /* sigma_8 */
	      params_trans[4] = wmap_likely(wmap, finalstates+ic, err);
	      forwardError(*err, __LINE__,);
	      fprintf(stderr, "sigma_8 = %g\n", params_trans[4]);
	      /* Omegab */
	      params_trans[5] = finalstates[ic+5]/100.0/dsqr(finalstates[ic+1]);
	   }
	   break;

	 default :
	    *err = addError(wmap_unknown, "unknown parametrization", *err, __LINE__);
	    return;
      }

      for (j=0; j<Ntrans; j++) last[j] = finalstates[ic+j];

      print_step_lc(CHTRA, MC_AF_ACCEPT, 0.0, Ntrans, params_trans);

      /*
      sprintf(fname, "scalCls_%04d.dat", i);
      fprintf(stderr, "creating file %s\n", fname);
      rename("scalCls.dat", fname);
      */

   }


   free(params_trans);
			 free(last);
   free(wmap);

   fclose(CHTRA);
}
#undef eps

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


   read_config_file(&config, "config_mcmc", &err); exitOnError(err, stderr);
   if (config.data[0]!=CMB) {
      fprintf(stderr, "data set has to be CMB");
      exit(2);
   }

   sprintf(fname, "%s.%s", chainout, chainfin_suf);
   CHAIN = fopen(fname, "r");
   if (!CHAIN) {
      fprintf(stderr, "%s.%s not found", chainout, chainfin_suf);
      exit(mcmc_io);
   }

   finalstates = read_mkmc_chain(config.npar, config.nchain, CHAIN, stderr, &nfinal, &headpos, NULL, NULL, &err);
   exitOnError(err, stderr);
   fclose(CHAIN);


   transform_chain(config, finalstates, nfinal, &err); exitOnError(err, stderr);

   return 0;
}

	    /*
	      sigma_8  = dsqr(D_plus(model, 1.0, 0, err));  forwardError(*err, __LINE__,);
	      sigma_8 *= sigma_R_sqr(model, 8.0, err);      forwardError(*err, __LINE__,); 
	      params_trans[5] = sigma_8;
	    */

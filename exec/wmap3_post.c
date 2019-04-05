/* ============================================================ *
 * wmap3_post.c							*
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
#include "config.h"
#include "mvdens.h"

#include "param.h"
#include "mcmc.h"
#include "cosmo_mcmc.h"
#include "stdnames.h"

#include "wmap.h"
#include "lensing.h"
#include "C_wrappers.h"

#define NNZ 5
void transform_chain(config_mcmc config, const double *finalstates, int nfinal, error **err)
{
   int i, Ntrans, ic;
   wmap_like *wmap;
   double *params_trans, sigma_8;
   FILE *CHTRA;
   char fname[100];
   cosmo *model;
   double par_nz[NNZ] = {0.0, 6.0, 0.61, 8.13, 0.62};

   Ntrans = 7;

   wmap = (wmap_like*)config.data_extra[0];
   sprintf(fname, "%s.tra", chainout);
   CHTRA = fopen(fname, "w");
   testErrorRet(!CHTRA, mcmc_io, "Could not create chain output file", *err, __LINE__,);
   //write_header_lc(CHTRA, Ntrans);

   params_trans = vector(0, Ntrans, err);   forwardError(*err, __LINE__,);
   model = init_parameters(0.3, 0.7, -1.0, 0.0, 0.7, 0.04, 0.0, 0.0, 0.9, 1.0,
			   NNZ, par_nz, smith03,
			   eisenhu, growth_de, linder, ymmk, SCOCOU,
			   norm_s8, err);
   forwardError(*err, __LINE__,);

   for (i=0; i<nfinal; i++) {

      ic = i*config.npar;

      switch (config.param) {

	 case param_wCDM_flat_nphys :

	    /* omegab = Omegab h^2 */
	    params_trans[0] = finalstates[ic+0]*finalstates[ic+6]*finalstates[ic+6];
	    /* Omegam = Omegab + Omegac */
	    params_trans[1] = finalstates[ic+0]+finalstates[ic+1];
	    params_trans[2] = finalstates[ic+2];
	    params_trans[3] = finalstates[ic+3];
	    params_trans[4] = finalstates[ic+4];
	    params_trans[6] = finalstates[ic+6];

	    updateParameters(model,
			     params_trans[1],
			     1.0-params_trans[1],
			     params_trans[3]  ,
			     model->w1_de  ,
			     params_trans[6],
			     finalstates[ic],
			     model->Omega_nu_mass,
			     model->Neff_nu_mass,
			     exp(finalstates[ic+5])*1.0e-10,        /* As */
			     params_trans[4],
			     model->Nnz,
			     model->par_nz ,
			     model->nonlinear ,
			     model->transfer  ,
			     model->growth    ,
			     model->de_param  ,
			     model->nofz      ,
			     model->bispmode  ,
			     norm_as,
			     err);
	    forwardError(*err, __LINE__,);

	    /* sigma_8 */
	    /*
	    sigma_8  = dsqr(D_plus(model, 1.0, 0, err));  forwardError(*err, __LINE__,);
	    sigma_8 *= sigma_R_sqr(model, 8.0, err);      forwardError(*err, __LINE__,); 
	    params_trans[5] = sigma_8;
	    */
	    sigma_8 = mysigma_8(finalstates[ic], finalstates[ic+1], 1.0-finalstates[ic]-finalstates[ic]+1,
				0.0, finalstates[ic+6]*100, finalstates[ic+5], finalstates[ic+4],
				finalstates[ic+3], finalstates[ic+2]);
	    dump_param(model, stderr);
	    fprintf(stderr, "sigma_8 = %g\n", sigma_8);
	    exit(0);
	    break;

	 default :
	    *err = addError(wmap_unknown, "unknown parametrization", *err, __LINE__);
	    return;
      }

      print_step_lc(CHTRA, MC_AF_ACCEPT, 0.0, Ntrans, params_trans);
   }


   free_vector(params_trans);
   free(wmap);

   fclose(CHTRA);
}
#undef NNZ

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

   finalstates = read_mkmc_chain(config.npar, config.nchain, CHAIN, stderr, &nfinal, &headpos, &err);
   exitOnError(err, stderr);
   fclose(CHAIN);


   transform_chain(config, finalstates, nfinal, &err); exitOnError(err, stderr);

   return 0;
}


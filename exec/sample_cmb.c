/* ============================================================ *
 * Martin Kilbinger 2008					*
 * Reads a sample (MCM chain or PMC (re)sample and prints the   *
 * corresponding Cl and Pdelta for each parameter vector.	*
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

#include "wmap.h"
#include "C_wrappers.h"

int same_parameter(const double *X1, const double *X2, int npar)
{
   int i;

   for (i=0; i<npar; i++) {
      if (fabs(X1[i]-X2[i])>EPSILON1) return 0;
   }
   return 1;
}

int main(int argc, char *argv[])
{
   config_pmc config;
   error *myerr = NULL, **err;
   parabox *pb;
   pmc_simu *psim;
   int i, prev_i;
   FILE *F;
   wmap_state *wmapstate;
   common_like *like;
   double logL, prev_logL;
   char str[128], str2[128];

   err = &myerr;

   if (argc!=2) {
      fprintf(stderr, "Usage: sample_cmb sample-file\n");
      return 2;
   }

   /* Reading header of config file */
   read_config_pmc_file(&config, "config_pmc", NULL, NULL, 1, err);
   quitOnError(*err, __LINE__, stderr);

   /* Initialise parameter box */
   pb = parabox_from_config(config.base.npar, config.base.min, config.base.max, err);
   quitOnError(*err, __LINE__, stderr);

   /* Initialise PMC simulation */
   F        = fopen_err(argv[1], "r", err);
   quitOnError(*err, __LINE__, stderr);
   /* NEW: with clip weights */
   psim = pmc_simu_from_file(F, config.nsamples, config.base.npar, config.base.n_ded, NULL, config.nclipw, err);
   quitOnError(*err, __LINE__, stderr);
   fclose(F);


   /* Initialise CAMB and WMAP */

   /* Initialise a new wmap state variable, do not take the one from the config file */
   wmapstate = malloc_err(sizeof(wmap_state), err);
   quitOnError(*err, __LINE__, stderr);
   /* MKDEBUG TODO: set paths, lmax */
   sprintf(wmapstate->scamb_path, "%s", "scamb");
   sprintf(wmapstate->data_path, "%s", "data");

   like = extra_init_common(CMB, config.base.par, config.base.npar, err);
   quitOnError(*err, __LINE__, stderr);
   init_CMB(like, err);
   quitOnError(*err, __LINE__, stderr);

   prev_logL = -1.0;
   prev_i    = 0;

   for (i=0; i<psim->nsamples; i++) {
      if (i>0 && same_parameter(psim->X+i*config.base.npar, psim->X+prev_i*config.base.npar, config.base.npar)==1) {
         logL = prev_logL;
         printf("Same parameter\n");

         sprintf(str, "scalCls_%05d.dat", i);
	 sprintf(str2, "scalCls_%05d.dat", prev_i);
         symlink(str2, str);

         sprintf(str, "lensedCls_%05d.dat", i);
	 sprintf(str2, "lensedCls_%05d.dat", prev_i);
         symlink(str2, str);

         sprintf(str, "Pdelta_%05d.dat", i);
	 sprintf(str2, "Pdelta_%05d.dat", prev_i);
         symlink(str2, str);
      } else {
         logL = likeli_CMB(like, psim->X+i*config.base.npar, err);
         quitOnError(*err, __LINE__, stderr);
         printf("%d logL = %g\n", i, logL);

         sprintf(str, "scalCls_%05d.dat", i);
         rename("scalCls.dat", str);

         sprintf(str, "lensedCls_%05d.dat", i);
         rename("lensedCls.dat", str);

         sprintf(str, "Pdelta_%05d.dat", i);
         rename("Pdelta.dat", str);

         prev_logL = logL;
         prev_i    = i;
      }
   }
   
   pmc_simu_free(&psim);
   free_parabox(&pb);
   free(wmapstate);

   return 0;
}

/* ============================================================ *
 * add_deduced_halomodel.c					*
 * Martin Kilbinger 2009					*
 *								*
 * ============================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>


#include "pmclib/pmc.h"
#include "pmctools/errorlist.h"
#include "pmctools/io.h"
#include "pmctools/maths.h"

#include "nicaea/halomodel.h"

#include "param.h"
#include "halo.h"
#include "types.h"
#include "exec_helper.h"

double calculate_deduced_halomodel(cosmo_hm *hm, par_t par, double zm, error **err)
{
   double p_ded, z, a;

   p_ded = 0.0;

   if (zm < 0) {
      z = zmean(hm->redshift, 0, err);
      forwardError(*err, __LINE__, -1.0);
   } else {
      z = zm;
   }
   a = 1.0/(1.0+z);

   switch (par) {
     
      case p_Mhalo_av :
        p_ded = av_halo_mass(hm, a, err);forwardError(*err, __LINE__, -1.0);
        break;
      case p_log10Mhalo_av :
        p_ded = av_halo_mass(hm, a, err);forwardError(*err, __LINE__, -1.0);
        p_ded = log10(p_ded);
        break;
      case p_bgal_av :
        p_ded = av_gal_bias(hm, a, err); forwardError(*err, __LINE__, -1.0);
        break;
      case p_Ngal_av :
        p_ded = Ngal_mean(hm, a, err);   forwardError(*err, __LINE__, -1.0);
        break;
      case p_fr_sat :
        p_ded = av_frsat(hm, a, err);   forwardError(*err, __LINE__, -1.0);
        break;
      case p_ngal_den :
        p_ded = ngal_den(hm, a, logMmax, hm->log10Mstar_min, hm->log10Mstar_max, err); forwardError(*err, __LINE__, -1.0);
        break;
      case p_log10ngal_den :
        p_ded = ngal_den(hm, a, logMmax, hm->log10Mstar_min, hm->log10Mstar_max, err);  forwardError(*err, __LINE__, -1.0);
        p_ded = log10(p_ded);
        break;
      default :
        *err = addErrorVA(hm_par, "Parameter %d (%s) is not a valid deduced parameter",
		          *err, __LINE__, par, spar_t(par));
        return -1.0;
      }
   
   return p_ded;
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

   fprintf(stderr, "Usage: add_deduced_halomodel [OPTIONS] PSIM [PAR_1 [PAR_2 [...]]]\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: config_pmc)\n");
   fprintf(stderr, "  -o OUTNAME       Ouput pmcsim name (default: psim+ded)\n");
   fprintf(stderr, "  PSIM             pmc simulation file (pmcsim_iter)\n");
   fprintf(stderr, "  PAR_i            String for deduced parameter #i. If not given, deduced\n");
   fprintf(stderr, "                    parameters are read from the config file (default)\n");

   if (ex>=0) exit(ex);
}

#define NSTR 128
int main(int argc, char *argv[])
{
   FILE *F;
   config_pmc config;
   error *myerr = NULL, **err;
   pmc_simu *psim;
   int i, j, n_ded, c, quiet;
   double *X_ded, zm;
   cosmo_hm *hm, *hmbase;
   char *cname, *outname;
   extern char *optarg;
   extern int optind, optopt;
   par_t *par;

   err = &myerr;

   cname   = NULL;
   outname = NULL;
   quiet   = 0;

   while ((c = getopt(argc, argv, ":c:o:qh")) != -1) {

      switch (c) {
	 case 'c' :
	    cname = optarg;
	    break;
	 case 'o' :
	    outname = optarg;
	    break;
         case 'q' :
            quiet = 1;
            break;
	 case 'h' :
	    usage(0, NULL);
	 case ':' :
	    usage(2, "Argument for -%c option missing\n", optopt);
	 case '?' :
	    usage(3, "Unknown option -%c\n\n", optopt);
      }

   }

   if (cname==NULL) {
      cname = malloc_err(NSTR*sizeof(char), err);
      quitOnError(*err, __LINE__, stderr);
      strcpy(cname, "config_pmc");
   }

   read_config_pmc_file(&config, cname, NULL, NULL, 0, err);
   quitOnError(*err, __LINE__, stderr);

   if (argc - 1 - optind == 0) {

      /* Get parameters from config file */
      n_ded = config.base.n_ded;
      if (n_ded == 0) {
         usage(4, "Error: No deduced parameter found on command line and config file\n");
      }
      if (! quiet) fprintf(stderr, "Reading %d deduced parameter%s from config file\n", n_ded, n_ded==1 ? "": "s");
      par = config.base.par + config.base.npar;

   } else {

      /* Get parameters from command line */
      n_ded = argc - 1 - optind;
      if (! quiet) {
         fprintf(stderr, "Reading %d deduded parameter%s from command line", n_ded, n_ded==1 ? "": "s");
         if (config.base.n_ded > 0) fprintf(stderr, ", overriding config file");
         fprintf(stderr, "\n");
      }
      spar_to_par(&par, n_ded, (const char**)(argv+optind+1), err);
      quitOnError(*err, __LINE__, stderr);

   }

   F = fopen_err(argv[optind], "r", err);
   quitOnError(*err, __LINE__, stderr);

   /* Read pmc sample file with n_ded set to zero, to avoid error message. */
   psim = pmc_simu_from_file(F, config.nsamples*config.fsfinal, config.base.npar, 0,
			     NULL, config.nclipw, err);
   quitOnError(*err, __LINE__, stderr);
   fclose(F);

   X_ded = malloc_err(sizeof(double)*psim->nsamples*(psim->n_ded+n_ded), err);
   quitOnError(*err, __LINE__, stderr);

   /* Copy previous deduced parameters (if any) */
   for (i=0; i<psim->nsamples; i++) {
      for (j=0; j<psim->n_ded; j++) {
	 X_ded[i*(psim->n_ded+n_ded)+j] = psim->X_ded[i*psim->n_ded+j];
      }
   }

   /* Get halomodel from the config info */
   hmbase = NULL;
   for (i=0; i<config.base.ndata; i++) {
      if (config.base.data[i]==GalCorr) {
	 hmbase = ((halo_state*)(config.base.data_extra[i]->state))->model;
	 break;
      }
   }
   if (hmbase==NULL) {
      fprintf(stderr, "Erorr: GalCorr data type not found\n");
      exit(4);
   }
   zm = zmean(hmbase->redshift, 0, err);
   quitOnError(*err, __LINE__, stderr);

   hm = copy_parameters_hm(hmbase, err);
   quitOnError(*err, __LINE__, stderr);

   /* Calculate deduced parameters */
   for (i=0; i<psim->nsamples; i++) {

      /* Update halomodel with current sample parameters */
      fill_parameters_hm(hm, psim->X+i*psim->ndim, config.base.par, config.base.npar, err);
      quitOnError(*err, __LINE__, stderr);
      updateFrom_hm(hmbase, hm, err);
      quitOnError(*err, __LINE__, stderr);

      //dump_param_hm(hm, stderr, err); quitOnError(*err, __LINE__, stderr);

      for (j=0; j<n_ded; j++) {
         /* Note: replace zm with -1 if mean redshift changes between samples */
	 X_ded[i*(psim->n_ded+n_ded)+j] = calculate_deduced_halomodel(hm, par[j], zm, err);
	 quitOnError(*err, __LINE__, stderr);
      }
   }

   psim->n_ded += n_ded;
   psim->X_ded  = X_ded;

   if (outname==NULL) {
      outname = malloc_err(sizeof(char)*128, err);
      quitOnError(*err, __LINE__, stderr);
      sprintf(outname, "%s+ded", argv[optind]);
   }
   out_pmc_simu_cosmo_pmc(outname, psim, config.base.par, 1.0, err);
   quitOnError(*err, __LINE__, stderr);

   /* Write config file with deduced parameters */
   /*
   resize_err(config.base.spar, config.base.npar*sizeof(char)*50,
	      (config.base.npar+n_ded)*sizeof(char)*50, 1, err);
   config.base.n_ded += n_ded;
   F = fopen_err("config_pmc+ded", "w", err);  
   quitOnError(*err, __LINE__, stderr);
   write_config_pmc_file(F, config, err);
   quitOnError(*err, __LINE__, stderr);
   fclose(F);
   exit(0);
   */
   return 0;
}

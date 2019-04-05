/* ============================================================ *
 * add_deduced_cosebis.c					*
 * Martin Kilbinger 2013                          					*
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

#include "nicaea/lensing.h"

#include "lens.h"
#include "types.h"
#include "exec_helper.h"
#include "param.h"



void usage(int ex, const char* str, ...)
{
   va_list ap;

   va_start(ap, str);
   if (str!=NULL) {
      vfprintf(stderr, str, ap);
      fprintf(stderr, "\n");
   }
   va_end(ap);

   fprintf(stderr, "Usage: add_deduced_cosebis [OPTIONS] PSIM\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: config_pmc)\n");
   fprintf(stderr, "  -o OUTNAME       Ouput pmcsim name (default: psim+ded)\n");
   fprintf(stderr, "  -s               Orthogonal (uncorrelated) data vector on output\n");
   fprintf(stderr, "  PSIM             pmc simulation file (pmcsim_iter)\n");

   if (ex>=0) exit(ex);
}

#define NSTR 128
int main(int argc, char *argv[])
{
   FILE *F;
   config_pmc config;
   error *myerr = NULL, **err;
   pmc_simu *psim;
   int i, j, n_ded, c, quiet, orthogonal, nrot, k;
   double *X_ded;
   cosmo_lens *lens, *lensbase;
   lens_state *lensstate;
   char *cname, *outname;
   extern char *optarg;
   extern int optind, optopt;
   double *E_cos, *eigenvalue, **eigenvector, *y, sum_rms;
   //double chi2_o, chi2_t, *cinv;
   unsigned long int *sind;
   datcov *dc;

   err = &myerr;

   cname   = NULL;
   outname = NULL;
   orthogonal = 0;
   quiet   = 0;

   while ((c = getopt(argc, argv, ":c:o:sqh")) != -1) {

      switch (c) {
         case 'c' :
            cname = optarg;
            break;
         case 'o' :
            outname = optarg;
            break;
         case 's' :
            orthogonal = 1;
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

   /* Get lens model from the config info */
   lensbase = NULL;
   for (i=0; i<config.base.ndata; i++) {
      if (config.base.data[i]==Lensing) {
         lensstate = (lens_state*)(config.base.data_extra[i]->state);
         lensbase  = lensstate->model;
         break;
      }
   }
   if (lensbase==NULL) {
      fprintf(stderr, "Erorr: 'Lensing' data type not found\n");
      exit(4);
   }

   n_ded = lensstate->data->Ntheta;

   F = fopen_err(argv[optind], "r", err);
   quitOnError(*err, __LINE__, stderr);

   /* Read pmc sample file with n_ded set to zero, to avoid error message. */
   psim = pmc_simu_from_file(F, config.nsamples*config.fsfinal, config.base.npar, config.base.n_ded,
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

   
   lens = copy_parameters_lens_only(lensbase, err);
   quitOnError(*err, __LINE__, stderr);

   E_cos = malloc_err(sizeof(double) * n_ded, err);
   quitOnError(*err, __LINE__, stderr);

   y = malloc_err(sizeof(double) * n_ded, err);
   quitOnError(*err, __LINE__, stderr);

   if (orthogonal) {

      eigenvalue  = sm2_vector(1, n_ded, err);
      quitOnError(*err, __LINE__, stderr);
      eigenvector = sm2_matrix(1, n_ded, 1, n_ded, err);
      quitOnError(*err, __LINE__, stderr);

      /* Read covariance matrix again, at this point it is already in Cholesky-form */
      dc = init_datcov_for_cov_only(1, n_ded, err);
      quitOnError(*err, __LINE__, stderr);
      read_cov_tomo(dc, lensstate->covname_ptr[0], 0, err);  

      /* For chi^2 test */
      /*cinv = copy_matrix(dc->cov[0], n_ded, err);
      quitOnError(*err, __LINE__, stderr);
      sm2_inverse(cinv, n_ded, err);
      quitOnError(*err, __LINE__, stderr);*/

      jacobi_transform(dc->cov[0], n_ded, eigenvalue, eigenvector, &nrot, err);
      quitOnError(*err, __LINE__, stderr);

      /* Sort eigenvalues */
      sind  = malloc_err(sizeof(unsigned long int)*(n_ded+1), err);
      quitOnError(*err, __LINE__, stderr);
      indexx(n_ded, eigenvalue, sind, err);           
      quitOnError(*err, __LINE__, stderr);

      for (k=1,sum_rms=0.0; k<=n_ded; k++) {
         sum_rms += 1/sqrt(eigenvalue[sind[k]]);
      }

      /*
      printf("# Eigenvalues (increasing), rms, frac_rms:\n");
      for (k=1; k<=n_ded; k++) {
         printf("%d %g %g  %f\n", k, eigenvalue[sind[k]], sqrt(eigenvalue[sind[k]]), 1/sqrt(eigenvalue[sind[k]])/sum_rms);
      }
      printf("# Corresponding eigenvectors\n");
      for (k=1; k<=n_ded; k++) {
         printf("%d(%g)", k, eigenvalue[sind[k]]);
         for (j=1; j<=n_ded; j++) {
            printf(" % g", eigenvector[j][sind[k]]);
         }
         printf("\n");
      }
      */


      for (k=0; k<n_ded; k++) {
         /* Transform to orthogonal vectors: Transpose eigenvector matrix (eigenvectors are in columns) */
         y[k] = 0.0;
         for (j=0; j<n_ded; j++) {
            //printf("%g * %g", eigenvector[j+1][sind[k+1]], lensstate->data->data[j]);
            y[k] += eigenvector[j+1][sind[k+1]] * lensstate->data->data[j];
            //if (j!=n_ded-1) printf(" + ");
         }
         //printf(" = %g\n", y[k]);
      }

		if (!quiet) {
			printf("# Orthgonal modes corresponding to data file '%s'\n", lensstate->datname);
			printf("# k mode var COSEBIs\n");
			for (k=1; k<=n_ded; k++) {
				printf("%d %g %g  %g\n", k, y[k-1], sqrt(eigenvalue[sind[k]]), lensstate->data->data[k-1]);
			}
		}

   }

   /* Calculate deduced parameters */
   for (i=0; i<psim->nsamples; i++) {

      /* Update lensmodel with current sample parameters */
      fill_parameters_lens(lens, NULL, lensstate, psim->X+i*psim->ndim, config.base.par, config.base.npar, err);
      quitOnError(*err, __LINE__, stderr);
      updateFrom_lens(lensbase, lens, err);
      quitOnError(*err, __LINE__, stderr);

      //dump_param_lens(lens, stderr, 0, err); quitOnError(*err, __LINE__, stderr);

      /* Calculate COSEBIs */
      for (j=0; j<n_ded; j++) {
         E_cos[j] = E_cosebi(lens, j+1, lensstate->cosebi_info->th_min,
                             lensstate->cosebi_info->th_max, 0, 0, lensstate->cosebi_info->path, NULL, err);
         quitOnError(*err, __LINE__, stderr);
      }

      for (k=0; k<n_ded; k++) {
         /* Transform to orthogonal vectors */
         if (orthogonal) {
            y[k] = 0.0;
            for (j=0; j<n_ded; j++) {
               y[k] += eigenvector[j+1][sind[k+1]] * E_cos[j];
            }
         } else {
            y[k] =  E_cos[k];
         }
      }

      /* Copy to deduced parameter vector */
      for (j=0; j<n_ded; j++) {
         X_ded[i*(psim->n_ded+n_ded)+j] = E_cos[j];
      }

      /* Chi^2 test */
      /*
      for (k=0,chi2_o=0.0; k<n_ded; k++) {
         for (j=0; j<n_ded; j++) {
            chi2_o += E_cos[k] * cinv[k*n_ded + j] * E_cos[j];
         } 
      }
      for (k=0,chi2_t=0.0; k<n_ded; k++) {
         chi2_t += y[k] / eigenvalue[sind[k+1]] * y[k];
      } 
      printf("# chi^2_null (orig trans) = %g %g\n", chi2_o, chi2_t);
      */


   }

   free(E_cos);
   free(y);

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

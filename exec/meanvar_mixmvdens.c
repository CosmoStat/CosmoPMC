/* ============================================================ *
 * meanvar_mixmvdens.c						*
 * Martin Kilbinger 2008					*
 *								*
 * Calculates the mean and confidence intervals of a mix_mvdens *
 * distribution (Gaussian or Student-t mixture).		*
 * ============================================================ */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "pmclib/pmc.h"
#include "param.h"
#include "errorlist.h"
#include "io.h"
#include "math.h"
#include "exec_helper.h"


#define JMAX 40
double rtbis(double (*func)(double, const mix_mvdens *, const double *, int, int, const double *, int),
	     double x1, double x2, double xacc, const mix_mvdens *mixmvd, const double *mu,
	     int i, int s, const double *sigma, int sign, int verbose)
{
   int j;
   double dx,f,fmid,xmid,rtb;

   f    = (*func)(x1, mixmvd, mu, i, s, sigma, sign);
   fmid = (*func)(x2, mixmvd, mu, i, s, sigma, sign);
   //fprintf(stderr, "rtbis: %f %f | %f %f\n", x1, x2, f, fmid);
   if (f*fmid >= 0.0) {
      if (verbose) fprintf(stderr, "Root not bracketed, f(%f)=%f, f(%f)=%f (sign=%d)\n", x1, f, x2, fmid, sign);
      return -10.0;
      //exit(3);
   }

   rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
   for (j=1;j<=JMAX;j++) {
      fmid = (*func)(xmid=rtb+(dx *= 0.5), mixmvd, mu, i, s, sigma, sign);
      //fprintf(stderr, "    %f %f\n", xmid, fmid);
      if (fmid <= 0.0) rtb=xmid;
      if (fabs(dx) < xacc || fmid == 0.0) return rtb;
   }

   fprintf(stderr, "Too many steps in rtbis\n");
   exit(2);
}
#undef JMAX

double area_one(double weight, double muminusmud, double y, double sigma)
{
   double denom, res;

   denom = sqrt(2.0)*sigma;
   res = weight*(erf((muminusmud+y)/denom) - erf(muminusmud/denom));

   //printf("area_one %g %g   %g  %g %g %g\n", muminusmud, denom, sigma, erf(muminusmud/denom), y, weight);

   return res;
}

double area(double y, const mix_mvdens *mixmvd, const double *mu, int i, int s, const double *sigma, int sign)
{
   int d;
   double res, x;

   for (d=0,res=0.0; d<mixmvd->ncomp; d++) {
      if (mixmvd->wght[d]>0) {
	 x    = sign*area_one(mixmvd->wght[d], mu[i]-mixmvd->comp[d]->mean[i], sign*y, fabs(sigma[i]));
	 res += x;
	 //fprintf(stderr, "d x = %d %g %g %g\n", d, x, mixmvd->comp[d]->mean[i], mu[i]);
      }  
   }

   return res - erf(s/sqrt(2));
}

double pdf(double x, const mix_mvdens *mixmvd, int i, const double *sigma)
{
   int d;
   double res, sig;

   for (d=0,res=0.0; d<mixmvd->ncomp; d++) {
      sig = fabs(sigma[i]);
      //sig = (mixmvd->comp[d]->std[i*mixmvd->ndim+i]);
      //res += 1.0/sqrt(pi)/fabs(sigma[i]) * exp(dsqr((x-mixmvd->comp[d]->mean[i])/sigma[i])/2.0);
      if (mixmvd->wght[d]>0) {
	 res += mixmvd->wght[d] * 1.0/sqrt(pi)/sig * exp(-dsqr((x-mixmvd->comp[d]->mean[i])/sig)/2.0);
      }
   }

   return res;
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

   fprintf(stderr, "Usage: meanvar_mixmvdens [OPTIONS] FILE\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: config_pmc)\n");
   fprintf(stderr, "  -v               Verbose, printing error messages\n");
   fprintf(stderr, "  FILE             Mixmvdens file (e.g. proposal_?)\n");

   if (ex>=0) exit(ex);
}

#define N 100
#define NSTR 128
int main(int argc, char *argv[])
{
   FILE *F;
   mix_mvdens *mixmvd;
   error *myerr = NULL, **err;
   int d, i, s, ndim, nrot, c, verbose;
   double *mean, *eigenvalue, **eigenvector, *conf, x;
   char *cname;
   config_pmc config;
   extern char *optarg;
   extern int optind, optopt;


   err = &myerr;


   /* Command line options */
   cname = NULL;
   verbose = 0;
   while ((c = getopt(argc, argv, ":c:vh")) != -1) {

      switch (c) {
	 case 'c' :
	    cname = optarg;
	    break;
	 case 'v' :
	    verbose = 1;
	    break;
	 case 'h' :
	    usage(0, NULL);
	 case ':' :
	    usage(1, "Argument for -%c option missing\n", optopt);
	 case '?' :
	    usage(2, "Unknown option -%c\n\n", optopt);
      }

   }
   if (cname==NULL) {
      cname = malloc_err(NSTR*sizeof(char), err);
      quitOnError(*err, __LINE__, stderr);
      strcpy(cname, "config_pmc");
   }
   if (argc-optind>1) usage(6, "Too many arguments");
   if (argc-optind<1) usage(7, "Mixmvdens file not given");


   F = fopen_err(argv[optind], "r", err);   quitOnError(*err, __LINE__, stderr);
   mixmvd = mix_mvdens_dwnp(F, err);        quitOnError(*err, __LINE__, stderr);
   fclose(F);

   ndim = mixmvd->ndim;

   mean = calloc_err(ndim, sizeof(double), err);
   quitOnError(*err, __LINE__, stderr);
   for (i=0; i<ndim; i++) {
      for (d=0; d<mixmvd->ncomp; d++) {
	 mean[i] += mixmvd->wght[d]*mixmvd->comp[d]->mean[i];
      }
   }

   eigenvalue  = sm2_vector(1, ndim, err);
   quitOnError(*err, __LINE__, stderr);
   eigenvector = sm2_matrix(1, ndim, 1, ndim, err);
   quitOnError(*err, __LINE__, stderr);
   //printf("Eigenvalues:\n");
   for (d=0; d<mixmvd->ncomp; d++) {
      if (mixmvd->wght[d]>0) {
	 jacobi_transform(mixmvd->comp[d]->std, ndim, eigenvalue, eigenvector, &nrot, err);
	 quitOnError(*err, __LINE__, stderr);
      }

      /*
      printf("Component %d ", d);
      for (i=0; i<ndim; i++) {
	 printf(" % .5f", eigenvalue[i+1]);
      }
      printf("\n");*/

   }

   //conf = sm2_matrix(0, ndim, 0, 5, err);
   conf = calloc_err(6*ndim, sizeof(double), err);
   quitOnError(*err, __LINE__, stderr);
   for (i=0; i<ndim; i++) {
      for (s=0; s<3; s++) {
	 conf[i*ndim+s]   = rtbis(area, -100, 100, 0.001, mixmvd, mean, i, s+1, eigenvalue+1, +1, verbose);
	 conf[i*ndim+s+3] = rtbis(area, -100, 100, 0.001, mixmvd, mean, i, s+1, eigenvalue+1, -1, verbose);
	 //printf("%d %d %f\n", i, s, conf[i][s]);
	 //exit(0);
      }
   }

   /* Reading header of config file */
   read_config_pmc_file(&config, cname, NULL, NULL, 1, err);
   quitOnError(*err, __LINE__, stderr);

   print_mean_sigma("", mean, conf, config.base.par, ndim, avg_mean, err);
   quitOnError(*err, __LINE__, stderr);
   for (i=0; i<ndim; i++) {
      //printf("%2d %10.5g", i, mean[i]);
      //for (s=0; s<3; s++) printf("   % 9.5f % 9.5f", conf[i][s], conf[i][s+3]);
      //printf("\n");
   }

   for (i=0; i<ndim; i++) {
      sprintf(cname, "prop_likeli1d_%d", i);
      F = fopen(cname, "w");
      for (x=config.base.min[i]; x<=config.base.max[i]; x+=(config.base.max[i]-config.base.min[i])/N) {
	 fprintf(F, "% g % g\n", x, pdf(x, mixmvd, i, eigenvalue+1));
      }
      fclose(F);
   }


   free(mean);
   free(conf);
   mix_mvdens_free(&mixmvd);
   sm2_free_vector(eigenvalue, 1, ndim);
   sm2_free_matrix(eigenvector, 1, ndim, 1, ndim);
   //sm2_free_matrix(conf, 0, ndim, 0, 2);

   return 0;
}

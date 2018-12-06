/* ============================================================ *
 * cl_one_sided.c						*
 * Martin Kilbinger 2009.					*
 * One-sided confidence levels.					*
 * ============================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "pmclib/pmc.h"
#include "param.h"
#include "errorlist.h"
#include "pmclib/tools.h"
#include "exec_helper.h"


void usage(int ex, const char* str, ...)
{
   va_list ap;

   va_start(ap, str);
   if (str!=NULL) {
      vfprintf(stderr, str, ap);
      fprintf(stderr, "\n");
   }
   va_end(ap);

   fprintf(stderr, "Usage: cl_one_sided [OPTIONS] sample\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: config_pmc)\n");
   fprintf(stderr, "  -i INDEX         Parameter index\n");
   fprintf(stderr, "  -d DIR           Direction (DIR=+1,-1)\n");
   fprintf(stderr, "  -v VALUE         Starting value\n");
   fprintf(stderr, "  -w WHICH         WHICH=0: 68%%,95%%,99.7%% c.l. (default)\n");
   fprintf(stderr, "                   WHICH=1: 68%%,90%%,95%% c.l.\n");
   fprintf(stderr, "  sample           PMC sample file\n");
   fprintf(stderr, "  The options -i INDEX, -d DIR and -v VALUE are required\n");

   if (ex>=0) exit(ex);
}

#define NSTR 128
int main(int argc, char *argv[])
{
   const double conf_123[3]  = {conf_68, conf_95, conf_99};
   const double conf_11h2[3] = {conf_68, conf_90, conf_95};

   pmc_simu *psim;
   error *myerr = NULL, **err;
   FILE *PMCSIM;
   char *cname;
   config_pmc config;
   double val, cl[3];
   const double *conf_xxx;
   int i, j, dir, w, c;
   extern char *optarg;
   extern int optind, optopt;


   err = &myerr;


   /* Command line options */
   cname = NULL;
   i     = -1;
   val   = -100.0;
   dir   = 0;
   w     = 0;
   while ((c = getopt(argc, argv, ":c:i:v:d:w:h")) != -1) {

      switch (c) {
	 case 'c' :
	    cname = optarg;
	    break;
	 case 'i' :
	    i = atoi(optarg);
	    break;
	 case 'v' :
	    val = atof(optarg);
	    break;
	 case 'd' :
	    dir = atoi(optarg);
	    break;
	 case 'w' :
	    w = atoi(optarg);
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

   printf("val = %g\n", val);

   if (i==-1)              usage(3, "Option -i INDEX missing");
   if (val<-99.0)          usage(4, "Option -v VALUE missing");
   if (dir!=+1 && dir!=-1) usage(5, "direction (-d DIR) has to be +1 or -1");
   if (argc-optind>1)      usage(6, "Too many arguments");
   if (argc-optind<1)      usage(7, "Sample file not given");
   if (w!=0 && w!=1)       usage(8, "Wrong argument WHICH (-w)");

   /* Read config file (without initialising proposal) */
   read_config_pmc_file(&config, cname, NULL, NULL, 1, err);
   quitOnError(*err, __LINE__, stderr);

   PMCSIM  = fopen_err(argv[optind], "r", err);    quitOnError(*err, __LINE__, stderr);
   psim    = pmc_simu_from_file(PMCSIM, config.nsamples*config.fsfinal, config.base.npar, config.base.n_ded, NULL,
				config.nclipw, err);
   quitOnError(*err, __LINE__, stderr);

   if (w==0) {
      conf_xxx = conf_123;
   } else {
      conf_xxx = conf_11h2;
   }
   confidence_level(psim->X, psim->weights, psim->flg, psim->nsamples, psim->ndim, i, val, dir, conf_xxx, cl, err);
   quitOnError(*err, __LINE__, stderr);
   
   printf("%2d", i);
   for (j=0; j<3; j++) {
      printf(" %g", cl[j]);
   }
   printf("\n");

   for (j=0; j<3; j++) {
      printf(" %s %c %.4g at %g%%\n", spar_t(config.base.par[i]),
	     dir==+1?'<':'>', val+dir*cl[j], conf_xxx[j]*100);
   }
 
   return 0;
}



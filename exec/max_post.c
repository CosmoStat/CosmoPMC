/* ============================================================ *
 * max_post.c							*
 * Martin Kilbinger 2010					*
 * ============================================================ */

#include "max_post.h"


void get_max_method(max_method_t *max_method, const char *optarg, error **err)
{
   int c;

   STRING2ENUM(*max_method, optarg, max_method_t, smax_method_t, c, Nmax_method_t, err);
   forwardError(*err, __LINE__,);
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

   fprintf(stderr, "Usage: max_post [OPTIONS]\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: config_max)\n");
   fprintf(stderr, "  -m [c|a|n]       Maximum-search method: 'a' (amoeba, default), 'c' (cg),\n");
   fprintf(stderr, "                    'n' (none; print posterior for fiducial parameter and exit)\n");
   fprintf(stderr, "  -t               Test maximum at the end\n");
   fprintf(stderr, "  -s SEED          Use SEED for random number generator. If SEED=-1 (default)\n");
   fprintf(stderr, "                    the current time is used as seed.\n");
   fprintf(stderr, "  -p               Prints the maximum-posterior model to the file 'model_maxlog'\n");
   fprintf(stderr, "  -q               Quiet mode\n");
   fprintf(stderr, "  -o FILEOUT       Output file [default:maxlogP]\n");
   fprintf(stderr, "  -h               This message\n");

   if (ex>=0) exit(ex);
}

double get_prior(config_base config, error **err)
{
   double chi2;
   int i;
   special_t special;

   chi2 = 0.0;
   for (i=0; i<config.ndata; i++) {
      chi2   = config.logpr_default;
      //printf("1 chi2 = %g\n", chi2);
      special = get_special(config.data[i], config.data_extra[i], err);
      forwardError(*err, __LINE__, 0.0);
      chi2  += prior_log_pdf_special(special, config.par, config.min, config.max, config.npar, err);
      //printf("2 chi2 = %g\n", chi2);
      forwardError(*err, __LINE__, 0.0);
   }

   return chi2;
}

/* ============================================================ *
 * Main program.						*
 * ============================================================ */

#define NSTR 128
int main(int argc, char *argv[])
{
   error *myerr = NULL, **err;
   config_max config;
   char *cname, fname[NSTR], MAXP_NAME[1000];
   time_t t_start;
   FILE *FLOG, *MAXP;
   int c, do_test_maximum, npar, i, j, nfunk, quiet, seed, want_model;
   extern char *optarg;
   extern int optind, optopt;
   parabox *pb;
   double maxP, *pstate, **simplex_start, lambda, *y, chi2;
   gsl_rng *rng;
   max_method_t max_method;


   err = &myerr;


   /* Command line options */
   cname           = NULL;
   do_test_maximum = 0;
   want_model      = 0;
   quiet           = 0;
   max_method      = max_amoeba;
   seed            = -1;
   strcpy(MAXP_NAME,   "maxlogP");

   while ((c = getopt(argc, argv, ":c:m:pts:o:qh")) != -1) {

     switch (c) {
     case 'c' :
       cname = optarg;
       break;
     case 'm' :
       get_max_method(&max_method, optarg, err);
       quitOnError(*err, __LINE__, stderr);
       break;
     case 't' :
       do_test_maximum = 1;
       break;
     case 'p' :
       want_model = 1;
       break;
     case 's' :
       seed = atoi(optarg);
       break;
     case 'q' :
       quiet = 1;
       break;
      case 'o' :
       strcpy(MAXP_NAME, optarg );
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
      strcpy(cname, "config_max");
   }


   fprintf(stderr, "max_post started\n");

   if(strcmp(MAXP_NAME, "maxlogP")){
       sprintf(fname, "%s_max_post.%s", log_name, MAXP_NAME);
     }else{
       sprintf(fname, "%s_max_post", log_name);
     }
   FLOG = fopen_err(fname, "w", err);  quitOnError(*err, __LINE__, stderr);
   fprintf(FLOG, "%s v%s compiled on %s %s\n", __FILE__, COSMO_PMC_VERSION, __DATE__, __TIME__);
      t_start = start_time(FLOG);


   fprintf(stderr, "Config file: %s\n", cname);
   fprintf(FLOG, "Config file: %s\n", cname);
   read_config_max_file(&config, cname, err);      quitOnError(*err, __LINE__, stderr);
   out_config_max(FLOG, config, err);              quitOnError(*err, __LINE__, stderr);
   fflush(FLOG);

   npar = config.base.npar;

   pb = parabox_from_config(npar, config.base.min, config.base.max, err);
   quitOnError(*err, __LINE__, stderr);

   pstate = malloc_err(npar*sizeof(double), err);
   quitOnError(*err, __LINE__, stderr);

   rng = init_random(seed, FLOG);

   set_start_parameter(pstate, config.start, config.fid, config.base, rng, FLOG, err);
   quitOnError(*err, __LINE__, stderr);
   
   maxP = MC_LOGBIG;

   switch (max_method) {

      case max_cg :

	 /* Conjugate-gradient method */
	 fprintf(stderr, "Looking for maximum posterior parameter using the ");
	 fprintf(stderr, "cg (conjuage-gradient) method\n");
	 maxP = max_logP(&config.base, pstate, pb, config.tolerance, quiet, err);
	 quitOnError(*err, __LINE__, stderr);
	 break;

      case max_amoeba :

	 /* Simplex method (amoeba) */
	 fprintf(stderr, "Looking for maximum posterior parameter using the ");
	 fprintf(stderr, "amoeba (simplex) method\n");
	 simplex_start = sm2_matrix(0, npar, 0, npar-1, err);
	 quitOnError(*err, __LINE__, stderr);
	 /* First simplex vector = starting vector */
	 for (i=0; i<npar; i++) {
	    simplex_start[0][i] = pstate[i];
	 }
	 /* Further: j-th point: add fraction of boxsize in (j-1) direction */
	 for (j=1; j<=npar; j++) {
	    for (i=0; i<npar; i++) {
	       simplex_start[j][i] = simplex_start[0][i];
	    }
	    lambda = (config.base.max[j-1]-config.base.min[j-1])/5.0;
	    simplex_start[j][j-1] += lambda;
	 }
	 y = amoeba_starting_values(simplex_start, npar, &config.base, pb, quiet, err);
	 quitOnError(*err, __LINE__, stderr);
	 amoeba(simplex_start, y, npar, config.tolerance, &nfunk, &config.base, pb, quiet, err);
	 quitOnError(*err, __LINE__, stderr);

	 /* Find vertex of final simplex with highest logP */
	 for (j=0,c=-1; j<=npar; j++) {
	    if (y[j]<maxP) {
	       maxP = y[j];
	       c    = j;
	    }
	 }
	 if (c<0) {
	    fprintf(stderr, "Error: Maximum-posterior point not found in final simplex.\n");
	    fprintf(stderr, "If no other obvious error occured (see screen printout above):\n");
	    fprintf(stderr, "Try using a different (fiduical) starting point (key 'fid' in 'config_max')\n");
	    return 4;
	 }
	 
	 /* Copy best parameter vector */
	 for (i=0; i<npar; i++) {
	    pstate[i] = simplex_start[c][i];
	 }
	 break;

      case max_none :
	
	
	 maxP = post_cg(&config.base, pstate, pb, quiet, err);
	 quitOnError(*err, __LINE__, stderr);
	 break;

      default :
	 usage(3, "Wrong max_method type");
   }

   /* "-maxP" = logP */
   chi2 = -2.0 * (-maxP - get_prior(config.base, err));
   quitOnError(*err, __LINE__, stderr);

   //post_cg(&config.base, pstate, pb, err); quitOnError(*err, __LINE__, stderr);

   fprintf(FLOG, "max log P = %g ftol = %g p = ", -maxP, config.tolerance);
   print_parameter(FLOG, config.base.npar, pstate); fflush(FLOG);
   fprintf(FLOG, "min chi^2 = %g\n", chi2);

   fprintf(stderr, "max log P = %g ftol = %g p = ", -maxP, config.tolerance);
   print_parameter(stderr, config.base.npar, pstate);
   fprintf(stderr, "min chi^2 = %g\n", chi2);

   MAXP = fopen_err(MAXP_NAME, "w", err);
   quitOnError(*err, __LINE__, stderr);
   fprintf(MAXP, "max log P = %g ftol = %g p = ", -maxP, config.tolerance);
   print_parameter(MAXP, config.base.npar, pstate);
   fclose(MAXP);

   if (want_model) {
      write_mean(pstate, "model_maxlog", config.base.ndata, (void**)config.base.data_extra, err);
      quitOnError(*err, __LINE__, stderr);
   }

   if (do_test_maximum) {
      test_maximum(pstate, config.base, fh, -maxP, err);
      quitOnError(*err, __LINE__, stderr);
   }

   /* Clean up */
   gsl_rng_free(rng);
   free(pstate);
   free_parabox(&pb);
   if (max_method==max_amoeba) sm2_free_matrix(simplex_start, 0, npar, 0, npar-1);

   fprintf(FLOG, "max_post ");
   end_time(t_start, FLOG);
   fclose(FLOG);

   fprintf(stderr, "max_post finished\n");

   return 0;
}

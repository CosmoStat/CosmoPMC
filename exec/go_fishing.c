/* ============================================================ *
 * go_fishing.c							*
 * Martin Kilbinger 2010					*
 * ============================================================ */

#include "go_fishing.h"

/* ============================================================ *
 * Returns the fisher matrix element (ab[0], ab[1]).		*
 * ============================================================ */

double fisher_element_adaptive(config_base *config, const double *pos, int ab[], int quiet, error **err)
{
   double *param;
   double h, f_ab, errn;
   int i;

   param = sm2_vector(1, config->npar, err); forwardError(*err, __LINE__, 0.0);
   for (i=0; i<config->npar; i++) {
      param[i] = pos[i];
   }

   h = 0.1;  // ad hoc

   /* f_ab = d^2[-logL]_{dab} */
   f_ab = -nd_dfridr2(posterior_log_pdf_common_void, ab[0], ab[1], param, h, h, (void*)config, &errn, err);
   forwardError(*err, __LINE__, 0.0);

   sm2_free_vector(param, 1, config->npar);

   if (! quiet) _DEBUGHERE_("fish(%d,%d) = %15.10f",ab[0], ab[1], -f_ab);

   return f_ab;
}

#define eps 1.0e-20
double fisher_element(config_base *config, const double *pos, int ab[], int quiet, error **err)
{
   double h[2], *param, f_ab, chi2[4];
   int diff[4][2] = {{+1,+1}, {+1,-1}, {-1,+1}, {-1,-1}};
   int i, j, k;

   param = sm2_vector(1, config->npar, err); forwardError(*err, __LINE__, 0);

   for (j=0; j<2; j++) {
      h[j] = fh*(config->max[ab[j]] - config->min[ab[j]]);                 /* fh times box size */
      testErrorRet(h[j]<eps, ce_infnan, "h too small", *err, __LINE__, 0);
   }

   /* Second derivative of the log-likelihood, Numerical Recipes (5.7.10) */
   for (j=0; j<4; j++) {

      /* Reset parameters */
      for (i=0; i<config->npar; i++) {
	 param[i] = pos[i];
      }

      /* Add small offsets */
      for (k=0; k<2; k++) param[ab[k]] += diff[j][k]*h[k];

      /* Symmetric on diagonal */
      if (j==2 && ab[0]==ab[1]) {
	 chi2[2] = chi2[1];
	 continue;
      }

      chi2[j] = posterior_log_pdf_common(config, param, err);

      //#ifdef DEBUG
      //fprintf(stderr, "fisher_element (%d,%d): log P = % .2f   par = ", ab[0], ab[1], -chi2[j]);
      //print_parameter(stderr, config->npar, param);
      //#endif

      forwardError(*err, __LINE__, 0);
   }

   /* F_ab = d[-logL]_{,ab} ! */
   f_ab = -(chi2[0] - chi2[1] - chi2[2] + chi2[3])/(4.0*h[0]*h[1]);

   sm2_free_vector(param, 1, config->npar);

   if (! quiet) _DEBUGHERE_("fish(%d,%d) = %15.10f",ab[0], ab[1], -f_ab);

   return f_ab;
}
#undef eps

/* Obsolete. Instead, fisher_matrix_part is used (for parallel calculation). */
mvdens *fisher_matrix(config_base *config, const double *pos, int quiet, error **err)
{
   int ab[2];
   mvdens *fish;
   double res;

   fish = mvdens_alloc(config->npar, err);    forwardError(*err, __LINE__, 0);

   for (ab[0]=0; ab[0]<config->npar; ab[0]++) {
      for (ab[1]=0; ab[1]<=ab[0]; ab[1]++) {

	 if (config->initial==mcmcini_fisher_diag && ab[0]!=ab[1])
	   continue;    /* Only calculate diagonal */

	 /* Add Fisher matrices of all experiments */
	 res = fisher_element(config, pos, ab, quiet, err);
	 forwardError(*err, __LINE__, 0);
	 fish->std[ab[0]*config->npar+ab[1]] = res;

      }
   }

   /* Symmetrise */
   for (ab[1]=1; ab[1]<config->npar; ab[1]++) {
      for (ab[0]=0; ab[0]<ab[1]; ab[0]++) {
	 fish->std[ab[0]*config->npar+ab[1]] = fish->std[ab[1]*config->npar+ab[0]];
      }
   }

   return fish;
}

/* ============================================================== *
 * Calculate part of Fisher matrix, elements my_start...my_end-1. *
 * ============================================================== */
void fisher_matrix_part(config_base *config, const double *pos, double *fish_part,
			int my_start, int my_end, int adaptive, int quiet, error **err)
{
   int ab[2], count, index;

   fprintfDEBUG(stderr, "fisher_matrix_part (elements %d .. %d)\n", my_start, my_end-1);

   for (ab[0]=index=count=0; ab[0]<config->npar; ab[0]++) {
      for (ab[1]=0; ab[1]<=ab[0]; ab[1]++) {

	 if (!(config->initial==mcmcini_fisher_diag && ab[0]!=ab[1])) {  /* Only calculate diagonal */

	    if (!(index<my_start || index>=my_end)) {

	       if (adaptive==0) {
		  fish_part[count] = fisher_element(config, pos, ab, quiet, err);
	       } else {
		  fish_part[count] = fisher_element_adaptive(config, pos, ab, quiet, err);
	       }
	       forwardError(*err, __LINE__,);

	       count++;

	    }

	    index++;

	 }

      }
   }

}

/* ============================================================ *
 * Recieves Fisher matrix parts from nodes and assembles the    *
 * entire matrix.						*
 * ============================================================ */
mvdens *assemble_fisher(int npar, int nproc, initial_t initial, int average_nel, int extra_nel,
			const double *fish_part_master, int my_nel_master, error **err)
{
   mvdens *fish;
   int c, my_nel, a, b, my_start, my_end, index, count;
   double *fish_part;

   fish_part = malloc_err((average_nel+1)*sizeof(double), err);     forwardError(*err, __LINE__, NULL);
   fish      = mvdens_alloc(npar, err);                             forwardError(*err, __LINE__, NULL);

   for (c=0; c<nproc; c++) {

      if (c != 0) {

	 /* Receive data from cliens */

	 /* Receive size */
	 MPI_Recv(&my_nel, 1, MPI_INT, c, gf_tag_size, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	 fprintfDEBUG(stderr, "Master: recieved my_nel=%d from node #%d\n", my_nel, c);

	 if (my_nel==0) continue;

	 /* Receive data */
	 MPI_Recv(fish_part, my_nel, MPI_DOUBLE, c, gf_tag_fisher_part, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	 fprintfDEBUG(stderr, "Master: recieved vector %p of size %d from node #%d\n", fish_part, my_nel, c);

      } else {

	 /* Use master data */

	 my_nel    = my_nel_master;
	 for (a=0; a<my_nel; a++) {
	    fish_part[a] = fish_part_master[a];
	 }

      }

      /* Copy to master Fisher matrix */
      my_start = c*average_nel + MIN(c, extra_nel);
      my_end   = my_start + my_nel;
      for (a=index=count=0; a<npar; a++) {
	 for (b=0; b<=a; b++) {

	    //fprintfDEBUG(stderr, "Testing2 (%d %d) (%d %d)\n", a, b, count, index);

	    if (!(initial==mcmcini_fisher_diag && a!=b)) {    /* Only calculate diagonal */

	       if (!(index<my_start || index>=my_end)) {

		  fish->std[a*npar+b] = fish_part[count];
		  fprintfDEBUG(stderr, "Copying fisher[%d][%d] = %g (%d %d)\n", a, b, fish_part[count], count, index);
		  count++;

	       }

	       index++;

	    }

	 }
      }

   }

   return fish;
}

void symmetrize_fisher(mvdens *fish)
{
   int a, b;

   for (b=1; b<fish->ndim; b++) {
      for (a=0; a<b; a++) {
	 fish->std[a*fish->ndim+b] = fish->std[b*fish->ndim+a];
      }
   }
}

/* Set off-diagonal to zero, absolute value on diagonal */
void force_positive(mvdens *fish)
{
   int a, b;

   for (a=0; a<fish->ndim; a++) {

      fish->std[a*fish->ndim + a] = fabs(fish->std[a*fish->ndim + a]);

      for (b=a+1; b<fish->ndim; b++) {
         fish->std[a*fish->ndim + b] = fish->std[b*fish->ndim + a] = 0.0;
      }
   }
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

   fprintf(stderr, "Usage: go_fishing [OPTIONS]\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c CONFIG        Configuration file (default: config_fish)\n");
   fprintf(stderr, "  -a               Adaptive numerical differentiation (default: fixed difference)\n");
   fprintf(stderr, "  -f               Force positive Fisher matrix\n");
   fprintf(stderr, "  -q               Quiet mode\n");
   fprintf(stderr, "  -h               This message\n");
   fprintf(stderr, "Run in parallel on NP cpu's: 'mpirun -np NP go_fishing [OPTIONS]\n");

   if (ex>=0) exit(ex);
}

/* ============================================================ *
 * Main program.						*
 * ============================================================ */

#define NSTR 128
int main(int argc, char *argv[])
{
   error *myerr = NULL, **err;
   config_fish config;
   char *cname, fname[NSTR];
   time_t t_start;
   FILE *FLOG;
   int c, i, myid, nproc, nel, my_nel, my_start, my_end, average_nel, extra_nel, adaptive, quiet, ret, force_pos;
   extern char *optarg;
   extern int optind, optopt;
   mvdens *fish;
   double *fish_part;


   err = &myerr;


   /* Command line options */
   cname     = NULL;
   adaptive  = 0;
   force_pos = 0;
   quiet     = 0;
   while ((c = getopt(argc, argv, ":c:afqh")) != -1) {

      switch (c) {
	 case 'c' :
	    cname = optarg;
	    break;
	 case 'a' :
	    adaptive = 1;
	    break;
	 case 'f' :
	    force_pos = 1;
	    break;
	 case 'q' :
	    quiet = 1;
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
      strcpy(cname, "config_fish");
   }


   /* Initialise  MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);

   if (myid==0) {
      fprintf(stderr, "go_fishing started\n");

      /* Log file */
      sprintf(fname, "%s_fish", log_name);
      FLOG = fopen_err(fname, "w", err);  quitOnError(*err, __LINE__, stderr);
      fprintf(FLOG, "%s v%s compiled on %s %s\n", __FILE__, COSMO_PMC_VERSION, __DATE__, __TIME__);
      t_start = start_time(FLOG);
      fprintf(FLOG, "Number of nodes nproc = %d\n", nproc);
      if (! quiet) {
	 fprintf(stderr, "%s v%s started on %d node%s\n", __FILE__, COSMO_PMC_VERSION, nproc, nproc>1?"s":"");
      }
   } else {
      FLOG = NULL;
   }


   /* Config file */
   sleep(myid*0.001);
   read_config_fish_file(&config, cname, err);      quitOnError(*err, __LINE__, stderr);
   if (myid==0) {
      out_config_fish(FLOG, config, err);              quitOnError(*err, __LINE__, stderr);
      fflush(FLOG);
   }


   if (myid==0 && ! quiet) {
      fprintf(stderr, "Calculating %s Fisher matrix...\n", config.base.initial==mcmcini_fisher_diag? "diagonal" : "");
      fflush(stderr);
   }


   /* Number of matrix elements to calculate */
   if (config.base.initial==mcmcini_fisher) {
      nel = config.base.npar*(config.base.npar+1)/2.0;
   } else if (config.base.initial==mcmcini_fisher_diag) {
      nel = config.base.npar;
   }

   average_nel = nel/nproc;   /* Zero if nproc>nel */
   extra_nel   = nel%nproc;

   if (! quiet) fprintfDEBUG(stderr, "nel nproc av extra = %d %d %d %d\n", nel, nproc, average_nel, extra_nel);

   /* Each node: initialise its part of the matrix */
   my_nel    = average_nel;
   if (myid<extra_nel) my_nel ++;   /* Distribute extra elements to first extra_nel nodes */

   /* Allocate memory with maximum size */
   fish_part = malloc_err((average_nel+1)*sizeof(double), err);        quitOnError(*err, __LINE__, stderr);
   my_start  = myid*average_nel + MIN(myid, extra_nel);
   my_end    = my_start + my_nel;


   /* All: Send size and data (NEW: master does not send to himself, did not work on sap/macbook) */
   if (myid != 0) {
      MPI_Send(&my_nel, 1, MPI_INT, 0, gf_tag_size, MPI_COMM_WORLD);
      //if (! quiet) fprintf(stderr, "proc %d >> sent my_nel = %d\n", myid, my_nel);
   }

   if (my_nel>0) {
      /* All: calculate part of Fisher matrix */
      if (! quiet) fprintf(stderr, "proc %d >> calling fisher_matrix_par (elements %d .. %d)\n", myid, my_start, my_end-1);
      fisher_matrix_part(&config.base, config.fid, fish_part, my_start, my_end, adaptive, quiet, err);
      quitOnError(*err, __LINE__, stderr);

      if (myid != 0) {
	 MPI_Send(fish_part, my_nel, MPI_DOUBLE, 0, gf_tag_fisher_part, MPI_COMM_WORLD);
	 //if (! quiet) fprintf(stderr, "proc %d >> sent fish_part vector %p of size %d\n", myid, fish_part, my_nel);
	 free(fish_part);
      }
   }


   /* Master: Receive parts and copy into Fisher matrix */
   if (myid==0) {

      fish = assemble_fisher(config.base.npar, nproc, config.base.initial, average_nel, extra_nel, 
			     fish_part, my_nel, err);
      quitOnError(*err, __LINE__, stderr);

      symmetrize_fisher(fish);

      if (force_pos == 1) {
	 force_positive(fish);
      }

      /* Mean */
      for (i=0; i<config.base.npar; i++) fish->mean[i] = config.fid[i];

      mvdens_chdump(fisher_name, fish, err);       quitOnError(*err, __LINE__, stderr);

      /* Testing positiveness */
      mvdens_cholesky_decomp(fish, err);
      if (isError(*err)) {
	 purgeError(err);
	 fprintf(stderr, "Fisher matrix not positive. You might not have found\n"
		 "the maximum-likelihood point. Possible solutions:\n"
		 " 1. Decrease fractional finite difference parameter ('fh' in 'mkmax.h') or\n"
       "    use adaptive differences (option '-a')\n"
		 " 2. Choose a fiducial start point closer to maximum (key 'fid' in 'config_fish'\n"
		 " 3. Only calculate diagonal Fisher matrix (set 'sinitial' to 'Fisher_diag' in 'config_fish'\n");
	 ret = 1;
      } else {
	 ret = 0;
      }

      mvdens_free(&fish);

      fprintf(stderr, "go_fishing finished\n");
      fprintf(FLOG, "go_fishing ");
      end_time(t_start, FLOG);
      fclose(FLOG);
   }


   MPI_Finalize();

   return ret;
}


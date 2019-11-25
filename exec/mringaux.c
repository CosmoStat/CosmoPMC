#include "mring.h"
#include "mvdens.h"
#include "config.h"
#include "param.h"

#if SANS_MRING == 0
#include "mringwr.h"
#endif

#include <sys/times.h>



void read_config_and_set_mring(const char cname[], config_max *config,
			       mring_state **mringstate, int *N, error **err)
{
   common_like *like;

   fprintf(stderr, "Config file: %s\n", cname);
   read_config_max_file(config, cname, err);
   forwardError(*err, __LINE__,);

   like = config->base.data_extra[0];
   *mringstate = (mring_state*)(like->state);
   //mring_print(stderr, (void*)(*mringlike));

   *N = like->npar + (*mringstate)->K;
}

double *read_a(const char aname[], int K, int N, error **err)
{
   FILE *F;
   double *a;
   int i;

   a = calloc_err(N, sizeof(double), err);      exitOnError(*err, stderr);

   /* Read coefficients a */
   fprintf(stderr, "Reading %d coefficients from %s\n", N, aname);
   F = fopen_err(aname, "r", err);              exitOnError(*err, stderr);
   for (i=K; i<N; i++) {
      fscanf(F, "%lg", a+i);
      fprintf(stderr, "a[%d]=%g ", i, a[i]);
   }
   fprintf(stderr, "\n");
   fclose(F);

   return a;
}

void F_test()
{
  int n; double  A, B;
  n = 0;
  A=1;
  B=1;

  printf("%f\n", Ff(n,A,B));
  printf("%f\n", Gf(n,A,B));

  
}

void S_N_Map(const char cname[])
{
   error *myerr = NULL, **err;
   cosmo_lens *lens;
   config_max config;
   mring_state *mringstate=NULL;
   double map2, cov, theta[2], res;
   int N;

   err = &myerr;

   read_config_and_set_mring(cname, &config, &mringstate, &N, err);
   exitOnError(*err, stderr);

   lens = set_cosmological_parameters_to_default_lens(err);
   exitOnError(*err, stderr);
   dump_param_lens(lens, stderr, 1, err);
   exitOnError(*err, stderr);

   map2 = map2_poly(lens, mringstate->THETA_MAX, 0, 0, err);
   exitOnError(*err, stderr);

   theta[0] = theta[1] = mringstate->THETA_MAX;
   cov = cov_Map2(theta, mringstate->cov_xi->theta, mringstate->cov_xi->cov, mringstate->cov_xi->n, err);
   exitOnError(*err, stderr);

   res = map2/sqrt(cov);

   printf("S_N_Map2: S/N = %f, Map^2 = %g, sqrt(cov) = %g, cov = %g\n", res, map2, sqrt(cov), cov);

   printf("%g\n", res);
   free_parameters_lens(&lens);
}

void S_N_Z_comb_none(const char cname[])
{
   double resZ, covR, eta, rr, dummy, theta[2], my_theta_min[2];
   error *myerr = NULL, **err;
   cosmo_lens *lens;
   config_max config;
   mring_state *mringstate=NULL;
   int N=0;

   err = &myerr;

   read_config_and_set_mring(cname, &config, &mringstate, &N, err);
   exitOnError(*err, stderr);

   eta      = mringstate->THETA_MIN/mringstate->THETA_MAX;
   theta[0] = theta[1] = mringstate->THETA_MAX;

   fprintf(stderr, "N=%d, theta_min=%g, theta_max=%g,  eta=%g  \n",
	   N, mringstate->THETA_MIN, mringstate->THETA_MAX, eta );

   lens = set_cosmological_parameters_to_default_lens(err);
   exitOnError(*err, stderr);
   dump_param_lens(lens, stderr, 1, err);
   exitOnError(*err, stderr);


   /* Using Z+ */
   R_from_xi(lens, mringstate->THETA_MIN/mringstate->THETA_MAX, mringstate->THETA_MAX, &rr, &dummy, err);
   exitOnError(*err, stderr);

   my_theta_min[0] = my_theta_min[1] = mringstate->THETA_MIN;
   covR = cov_RR_Z(my_theta_min, theta, mringstate->cov_xi->theta,
		   mringstate->cov_xi->cov, mringstate->cov_xi->n, err);
   exitOnError(*err, stderr);

   resZ = rr/sqrt(covR);
   printf("cov_RR_Z:   S/N = %f, RR = %g, sqrt(covR) = %g\n", resZ, rr, sqrt(covR));

   /* Print again */
   printf("%f %g\n", resZ, rr);


   free_parameters_lens(&lens);
}

void S_N_sum(const char cname[])
{
   error *myerr = NULL, **err;
   int n, N;
   double rr, covR, res, c0;
   cosmo_lens *lens;
   config_max config;
   mring_state *mringstate=NULL;
   double dtheta, dlogtheta, dummy, theta[2], my_theta_min[2];
   err = &myerr;


   read_config_and_set_mring(cname, &config, &mringstate, &N, err);
   exitOnError(*err, stderr);

   lens = set_cosmological_parameters_to_default_lens(err);
   exitOnError(*err, stderr);
   dump_param_lens(lens, stderr, 1, err);
   exitOnError(*err, stderr);


   if (mringstate->bintype==0) { /* linear bins */
      dtheta = (mringstate->THETA_MAX-mringstate->theta_max_start)/((double)(mringstate->Ntheta)-1.0);
      dlogtheta = 0.0;
   } else { /* logarithmis bins */
      dtheta = 0.0;
      dlogtheta = (log(mringstate->THETA_MAX)-log(mringstate->theta_max_start))/((double)(mringstate->Ntheta)-1.0);
   }

   theta[0] = mringstate->theta_max_start;
   res = 0.0;
   fprintf(stderr, " # theta_max s/n \n");
   for (n=0; n<mringstate->Ntheta; n++) {

      R_from_xi(lens, mringstate->THETA_MIN/theta[0], theta[0], &rr, &dummy, err);
      exitOnError(*err, stderr);

      theta[1] = theta[0];
      my_theta_min[0] = my_theta_min[1] = mringstate->THETA_MIN;
      covR = cov_RR_Z(my_theta_min, theta, mringstate->cov_xi->theta,
		      mringstate->cov_xi->cov, mringstate->cov_xi->n, err);
      exitOnError(*err, stderr);
      
      c0 = fabs(rr)/sqrt(covR);
      fprintf(stderr, "%lg %.2f \n",theta[0]/arcmin, c0);
      res += c0;
      if (mringstate->bintype==0) theta[0] += dtheta;
      else theta[0] *= exp(dlogtheta);
   }
   fprintf(stderr, "\n");
   res = res/(double)mringstate->Ntheta;

   printf("S/N (sum) for Z+ = %g\n", res);
}

void S_N_T(const char cname[], const char aname[])
{
   error *myerr = NULL, **err;
   int N=0, i;
   double *a, A, B, realA, realB, rr, covR, res, theta[2], tp, tm, y, x, eta, my_theta_min[2], sum;
   cosmo_lens *lens;
   config_max config;
   mring_state *mringstate=NULL;
   FILE *F;

   err = &myerr;

   read_config_and_set_mring(cname, &config, &mringstate, &N, err);
   exitOnError(*err, stderr);

   a = read_a(aname, mringstate->K, N, err);    exitOnError(*err, stderr);

   eta = mringstate->THETA_MIN/mringstate->THETA_MAX;
   realA = (mringstate->THETA_MAX - mringstate->THETA_MIN)/2.0;
   realB = (mringstate->THETA_MAX + mringstate->THETA_MIN)/2.0;
   if (mringstate->combinef!=all_sc) {
      A = realA;
      B = realB;
   } else {
      /* General optimisation, independent of scales, set A,B<0 */
      A = B = -1.0;
   }

   fill_constraint(a, N, mringstate->K, A, B, mringstate->poly, 0, err);
   exitOnError(*err, stderr);

   /* Normalisation for S/N ratio not necessary */
   normalise(a, N, mringstate->poly, err);
   exitOnError(*err, stderr);

   lens = set_cosmological_parameters_to_default_lens(err);
   exitOnError(*err, stderr);
   dump_param_lens(lens, stderr, 1, err);
   exitOnError(*err, stderr);

   fprintf(stderr, "\n\n\n");
   for (i=0,sum=0.0; i<N; i++) {
     fprintf(stderr, "a[%d] = % .5g;", i, a[i]);
     if (i==0 && mringstate->poly==cheby) sum += a[i]*a[i];
     else sum += a[i]*a[i]/2.0;
   }
   sum *= pi;
   fprintf(stderr, "\nSum = %g\n", sum);

   rr = RR(lens, mringstate->THETA_MIN, mringstate->THETA_MAX, a, N, mringstate->poly, +1, err);
   exitOnError(*err, stderr);
   
   theta[0] = theta[1] = mringstate->THETA_MAX;
   my_theta_min[0] = my_theta_min[1] = mringstate->THETA_MIN;
   covR = cov_RR(my_theta_min, theta, a, N, mringstate->poly,
		 mringstate->cov_xi->theta, mringstate->cov_xi->cov, mringstate->cov_xi->n,
		 outer, -1, -1, cov_xi_fac, err);
   exitOnError(*err, stderr);

   res = fabs(rr)/sqrt(covR);
   printf("S_N_T:     S/N = %f, RR = %g, sqrt(covR) = %g\n", res, rr, sqrt(covR));
   for (i=mringstate->K; i<N; i++) printf("%g ", a[i]);

   printf("\n%g\n", res);

   F = fopen("T+", "w");
   for (x=-1.0; x<=1.005; x+=0.01) {
      tp = Tp(x, a, N, mringstate->poly, err);
      exitOnError(*err, stderr);
      tm = Tm(x, a, N, mringstate->poly, realB/realA, err);
      exitOnError(*err, stderr);
      y = eta + (1.0-eta)/2.0*(1.0+x);
      fprintf(F, "%f % f % g % g\n", y, x, tp, tm);
   }
   fclose(F);

   free(a);
   free_parameters_lens(&lens);
}

void check_EB(const char cname[], const char aname[])
{
   error *myerr = NULL, **err;
   int N=0, i;
   double *a, A, B, realA, realB, rr_p, rr_m, eta, th_min;
   cosmo_lens *lens;
   config_max config;
   mring_state *mringstate=NULL;
   double x, tp, tm, y;
   //double c0, c1;
   FILE *F;

   err = &myerr;

 
   read_config_and_set_mring(cname, &config, &mringstate, &N, err);
   exitOnError(*err, stderr);

   eta = mringstate->THETA_MIN/mringstate->THETA_MAX;

   a = read_a(aname, mringstate->K, N, err);    exitOnError(*err, stderr);

   /* K constraints, N total coefficients (including constraints) */
   if (mringstate->THETA_MIN>=0) {
      realA = (mringstate->THETA_MAX - mringstate->THETA_MIN)/2.0;
      realB = (mringstate->THETA_MAX + mringstate->THETA_MIN)/2.0;
   } else {
      realA = (mringstate->THETA_MAX - mringstate->THETA_MAX*mringstate->eta)/2.0;
      realB = (mringstate->THETA_MAX + mringstate->THETA_MAX*mringstate->eta)/2.0;
   }
   if (mringstate->combinef!=all_sc) {
      A = realA;
      B = realB;
      eta = mringstate->eta;
   } else {
      /* General optimisation, independent of scales, set A,B<0 */
      A = B = -1.0;
      eta = -1.0;
   }


   /*
     // TESTING alpha_explicit
   x = -0.9343;
   i = 0;
   a1 = alpha_nointerp(x, i, N, realA, realB, mringstate->poly, err);
   exitOnError(*err, stderr);
   a2 = alpha_explicit(x, i, realA, realB, mringstate->poly, err);
   exitOnError(*err, stderr);
   printf("%d %g %g\n", i, a1, a2);

   i = 2;
   a1 = alpha_nointerp(x, i, N, realA, realB, mringstate->poly, err);
   exitOnError(*err, stderr);
   a2 = alpha_explicit(x, i, realA, realB, mringstate->poly, err);
   exitOnError(*err, stderr);
   printf("%d %g %g\n", i, a1, a2);
   //exit(0);
   */

   /* Use A, B only for fill_constraint! Otherwise realA, realB!! */

   lens = set_cosmological_parameters_to_default_lens(err);
   exitOnError(*err, stderr);

   //printf("^^^ %d %d %g %g %d  %f %f\n", N, mringstate->K, A/arcmin, B/arcmin, mringstate->poly,
   //  mringstate->THETA_MIN/arcmin, mringstate->THETA_MAX/arcmin);
   fill_constraint(a, N, mringstate->K, A, B, mringstate->poly, 0, err);
   exitOnError(*err, stderr);

   normalise(a, N, mringstate->poly, err);
   exitOnError(*err, stderr);


   for (i=0; i<N; i++) {
     fprintf(stderr, "a[%d] = % .10g \n", i, a[i]);
   }

   /*
   c1 = integrate_Tpm(1, a, realA, realB, N, mringstate->poly, +1, err);   exitOnError(*err, stderr);
   c0 = integrate_Tpm(0, a, realA, realB, N, mringstate->poly, +1, err);   exitOnError(*err, stderr);
   fprintf(stderr, "Integral constraints = %g %g", c0, c1);
   c1 = c0 = -1.0;
   //c1 = integrate_Tpm(1, a, realA, realB, N, mringstate->poly, -1, err);   exitOnError(*err, stderr);
   //c0 = integrate_Tpm(0, a, realA, realB, N, mringstate->poly, -1, err);   exitOnError(*err, stderr);
   fprintf(stderr, "   %g %g\n", c0, c1);

   for (i=0; i<mringstate->K; i++) {
      fprintf(stderr, "Constraint #%d = % g\n", i, check_constraint(a, i, N, mringstate->K, realA, realB,
								    mringstate->poly, err));
      exitOnError(*err, stderr);
   }
   */

   if (mringstate->THETA_MIN>0) th_min = mringstate->THETA_MIN;
   else th_min = mringstate->THETA_MAX*mringstate->eta;

   rr_p = RR(lens, th_min, mringstate->THETA_MAX, a, N, mringstate->poly, +1, err);
   exitOnError(*err, stderr);

   rr_m = RR(lens, th_min, mringstate->THETA_MAX, a, N, mringstate->poly, -1, err);
   exitOnError(*err, stderr);

   printf("%f %f  % g % g\n", th_min/arcmin, mringstate->THETA_MAX/arcmin,
	  0.5*(rr_p+rr_m), 0.5*(rr_p-rr_m));

   F = fopen("T+", "w");
   for (x=-1.0; x<=1.0; x+=0.01) {
      y = eta + (1.0-eta)/2.0*(1.0+x);
      tp = Tp(x, a, N, mringstate->poly, err);
      exitOnError(*err, stderr);
      tm = Tm(x, a, N, mringstate->poly, realB/realA, err);
      exitOnError(*err, stderr);
      fprintf(F, "%f %f % g % g\n", y, x, tp, tm);
   }
   fclose(F);

   free(a);
   free_parameters_lens(&lens);
}

void fisherT(const char cname[],  const char aname[], error **err)
{
   int n, N=0, diag_only;
   double *a, eta, A, B, res;
   double x, tp, tm, realA, realB, y;
   cosmo_lens *lens;
   config_max config;
   mring_state *mringstate=NULL;
   mvdens *fish;
   FILE *F;


   read_config_and_set_mring(cname, &config, &mringstate, &N, err);
   forwardError(*err, __LINE__,);

   if (!(mringstate->ltype==fisher_var || mringstate->ltype==fisher_covar)) {
      fprintf(stderr, "Not the correct ltype\n");
      exit(1);
   }

   a = read_a(aname, mringstate->K, N, err);      forwardError(*err, __LINE__,);

   /* K constraints, N total coefficients (including constraints) */
   if (mringstate->THETA_MIN>=0) {
      realA = (mringstate->THETA_MAX - mringstate->THETA_MIN)/2.0;
      realB = (mringstate->THETA_MAX + mringstate->THETA_MIN)/2.0;
   } else {
      realA = (mringstate->THETA_MAX - mringstate->THETA_MAX*mringstate->eta)/2.0;
      realB = (mringstate->THETA_MAX + mringstate->THETA_MAX*mringstate->eta)/2.0;
   }
   if (mringstate->combinef!=all_sc) {
      A = realA;
      B = realB;
      eta = mringstate->eta;
   } else {
      /* General optimisation, independent of scales, set A,B<0 */
      A = B = -1.0;
      eta = -1.0;
   }

   //printf("^^^ %d %d %g %g %d  %f %f\n", N, mringstate->K, A/arcmin, B/arcmin, mringstate->poly,
   //  mringstate->THETA_MIN/arcmin, mringstate->THETA_MAX/arcmin);
   fill_constraint(a, N, mringstate->K, A, B, mringstate->poly, 0, err);
   forwardError(*err, __LINE__,);

   normalise(a, N, mringstate->poly, err);
   forwardError(*err, __LINE__,);


   for (n=0; n<N; n++) {
      fprintf(stderr, "a[%d] = % .10f; ", n, a[n]);
   }
   fprintf(stderr, "\n");

   lens = set_cosmological_parameters_to_default_lens(err);
   forwardError(*err, __LINE__,);

   //printf("THETA_MIN eta A B theta_max_start = %g' %g  %g' %g'  %g'\n",
   //	  mringstate->THETA_MIN/arcmin, eta, A/arcmin, B/arcmin, mringstate->theta_max_start/arcmin);

   diag_only = mringstate->ltype==fisher_var ? 1 : 0;
   fish = fisher_matrix_mring(lens, mringstate->THETA_MIN, mringstate->theta_max_start, mringstate->THETA_MAX,
			      mringstate->Ntheta, mringstate->bintype, a, N, mringstate->poly, mringstate->par,
			      mringstate->npar, diag_only, mringstate->cov_xi, eta, NULL, err);
   forwardError(*err, __LINE__,);

   mvdens_chdump("fishT", fish, err);
   forwardError(*err, __LINE__,);

   if (mringstate->maxF==minusTrFm1 || mringstate->maxF==inverseTrFm1 ||
       mringstate->maxF==FoM || mringstate->maxF==FoMnocor) {
      sm2_inverse(fish->std, fish->ndim, err);
      forwardError(*err, __LINE__,);
   }

   mvdens_chdump("fishTinv", fish, err);
   forwardError(*err, __LINE__,);

   if (mringstate->maxF==FoM || mringstate->maxF==FoMnocor)  {
      /* Figure of merit = error ellipse area */
      res = figure_of_merit(fish, mringstate->maxF, err);
      forwardError(*err, __LINE__,);
   } else {
      /* Trace */
      for (n=0,res=0.0; n<fish->ndim; n++) res += fish->std[n*fish->ndim+n];
      if (mringstate->maxF==minusTrFm1) res = -res;
      else if (mringstate->maxF==inverseTrFm1) res = 1.0/res;
   }

   fflush(stdout);
   fprintf(stdout, "%g\n", res);


   F = fopen("T+", "w");
   for (x=-1.0; x<=1.001; x+=0.01) {
      tp = Tp(x, a, N, mringstate->poly,err);
      forwardError(*err, __LINE__,);
      tm = Tm(x, a, N, mringstate->poly, realB/realA, err);
      forwardError(*err, __LINE__,);
      y = eta + (1.0-eta)/2.0*(1.0+x);
      fprintf(F, "%f % f % g % g\n", y, x, tp, tm);
   }
   fclose(F);

   mvdens_free(&fish);
   free_parameters_lens(&lens);
   free(a);
}

/*****************************/
void S_N_random_a(const char cname[], const char aname[])
{
   error *myerr = NULL, **err;
   int i, N=0;
   double *a, A, B, realA, realB, rr, covR, res, theta[2], tp, tm, y, x, eta, my_theta_min[2];
   double ran, frac=0.1;
   struct tms time;
   clock_t cl;
   unsigned u;
  
   cosmo_lens *lens;
   config_max config;
   mring_state *mringstate=NULL;
   FILE *F;

   err = &myerr;

   read_config_and_set_mring(cname, &config, &mringstate, &N, err);
   exitOnError(*err, stderr);

   a = read_a(aname, mringstate->K, N, err);    exitOnError(*err, stderr);

   eta = mringstate->THETA_MIN/mringstate->THETA_MAX;
   realA = (mringstate->THETA_MAX - mringstate->THETA_MIN)/2.0;
   realB = (mringstate->THETA_MAX + mringstate->THETA_MIN)/2.0;
   if (mringstate->combinef!=all_sc) {
      A = realA;
      B = realB;
   } else {
      /* General optimisation, independent of scales, set A,B<0 */
      A = B = -1.0;
   }

   fill_constraint(a, N, mringstate->K, A, B, mringstate->poly, 0, err);
   exitOnError(*err, stderr);

   /* change all highest N-K a[n] randomly (1%) and  constraints to get the first K a[n] */

   /* random init */
   do {
      cl = times(&time);
      u = (unsigned)cl;
      u = (u%10)*(u%100)*(u%1000);
   } while (u==0);
   srand(u);

   for (i=mringstate->K;i<N;i++){
      ran=rand()/(double)RAND_MAX;
      a[i] = a[i] + (2*ran-1)*frac;
    }
 
   fill_constraint(a, N, mringstate->K, A, B, mringstate->poly, 0, err);
   exitOnError(*err, stderr);

   
   normalise(a, N, mringstate->poly, err);
   exitOnError(*err, stderr);
   /*
   for (i=0,sum=0.0; i<N; i++) {

      if (i==0 && mringstate->poly==cheby) sum += a[i]*a[i];
      else sum += a[i]*a[i]/2.0;
   }
   sum *= pi;
   fprintf(stderr, "\nSum = %g\n", sum);
   */

   
   lens = set_cosmological_parameters_to_default_lens(err);
   exitOnError(*err, stderr);
   dump_param_lens(lens, stderr, 1, err);
   exitOnError(*err, stderr);

   rr = RR(lens, mringstate->THETA_MIN, mringstate->THETA_MAX, a, N, mringstate->poly, +1, err);
   exitOnError(*err, stderr);
   
   theta[0] = theta[1] = mringstate->THETA_MAX;
   my_theta_min[0] = my_theta_min[1] = mringstate->THETA_MIN;
   covR = cov_RR(my_theta_min, theta, a, N, mringstate->poly,
		 mringstate->cov_xi->theta, mringstate->cov_xi->cov, mringstate->cov_xi->n,
		 outer, -1, -1, cov_xi_fac, err);
   exitOnError(*err, stderr);

   res = fabs(rr)/sqrt(covR);
   printf("S_N_T:     S/N = %f, RR = %g, sqrt(covR) = %g\n", res, rr, sqrt(covR));
   
   F = fopen("SN", "w");
   fprintf(F,"%g\n", res);
   fclose(F);
  

   F = fopen("T+", "w");
   for (x=-1.0; x<=1.005; x+=0.01) {
      tp = Tp(x, a, N, mringstate->poly, err);
      exitOnError(*err, stderr);
      tm = Tm(x, a, N, mringstate->poly, realB/realA, err);
      exitOnError(*err, stderr);
      y = eta + (1.0-eta)/2.0*(1.0+x);
      fprintf(F, "%f % f % g % g\n", y, x, tp, tm);
   }
   fclose(F);

   free(a);
   //   free_parameters_lens(&lens);
}

/****************************/





void extend_Tp_zero(const char cname[], const char aname[], double psi1)
{
   error *myerr = NULL, **err;
   config_max config;
   mring_state *mringstate=NULL;
   double *a0, *a1, A, B, x;
   int N=0, N1, n;
   FILE *F;

   err = &myerr;

   read_config_and_set_mring(cname, &config, &mringstate, &N, err);
   exitOnError(*err, stderr);

   if (psi1<mringstate->THETA_MAX) {
      fprintf(stderr, "Error: psi1 has to be larger than psi0\n");
      exit(1);
   }

   a0 = read_a(aname, mringstate->K, N, err);      exitOnError(*err, stderr);

   if (mringstate->combinef!=all_sc) {
      /* comb_none: mringstate->THETA_MIN is -1 */
      A = (mringstate->THETA_MAX - mringstate->THETA_MAX*mringstate->eta)/2.0;
      B = (mringstate->THETA_MAX + mringstate->THETA_MAX*mringstate->eta)/2.0;
   } else {
      /* all_sc */
      A = B = -1.0;
   }

   fill_constraint(a0, N, mringstate->K, A, B, mringstate->poly, 0, err);
   exitOnError(*err, stderr);

   //normalise(a0, N, mringstate->poly, err);
   //exitOnError(*err, stderr);

   for (n=0; n<N; n++) {
      fprintf(stderr, "a0[%d] = % .10f; ", n, a0[n]);
   }
   fprintf(stderr, "\n");

   N1 = N;
   a1 = decomp_Tp(N1, mringstate->eta, mringstate->poly, a0, mringstate->THETA_MAX, psi1, err);
   exitOnError(*err, stderr);

   for (n=0; n<N1; n++) {
      fprintf(stderr, "a1[%d] = % .10f; ", n, a1[n]);
   }
   fprintf(stderr, "\n");

   F = fopen_err("T+_extend", "w", err);
   exitOnError(*err, stderr);
   for (x=-1; x<=1; x+=0.01) {
      fprintf(F, "% f % f\n", x, Tp(x, a1, N1, mringstate->poly, err));
      exitOnError(*err, stderr);
   }
   fclose(F);

   F = fopen_err("list_a_extend", "w", err);
   exitOnError(*err, stderr);
   for (n=mringstate->K; n<N1; n++) {
      fprintf(F, "%.5f ", a1[n]);
   }
   fprintf(F, "\n");
   fclose(F);
}


void fisherZ(const char cname[], error **err)
{
   int N=0;
   double res, theta_min_ij[2], theta_max_ij[2], dtheta, dlogtheta, cov, eta;
   double theta_i, theta_j;
   int i, j, diag_only, n;
   cosmo_lens *lens;
   config_max config;
   mring_state *mringstate=NULL;
   mvdens *fish;
   FILE *F;


   /* N unused here */
   read_config_and_set_mring(cname, &config, &mringstate, &N, err);
   forwardError(*err, __LINE__,);

   lens = set_cosmological_parameters_to_default_lens(err);
   forwardError(*err, __LINE__,);

   if (mringstate->combinef!=all_sc) {
      eta = mringstate->eta;
   } else {
      /* General optimisation, independent of scales, set eta<0 */
      eta = -1.0;
   }

   /* Fisher matrix ... */
   diag_only = mringstate->ltype==fisher_var ? 1 : 0;
   printf("fisher_matrix_mring Z\n");
   fish = fisher_matrix_mring(lens, mringstate->THETA_MIN, mringstate->theta_max_start, mringstate->THETA_MAX,
			      mringstate->Ntheta, mringstate->bintype, NULL, N, mringstate->poly, mringstate->par,
			      mringstate->npar, diag_only, mringstate->cov_xi, eta, NULL, err);
   forwardError(*err, __LINE__,);
   printf("done\n");

   mvdens_chdump("fishZ", fish, err);
   forwardError(*err, __LINE__,);

   if (mringstate->maxF==minusTrFm1 || mringstate->maxF==inverseTrFm1 ||
       mringstate->maxF==FoM || mringstate->maxF==FoMnocor) {
      sm2_inverse(fish->std, fish->ndim, err);
      forwardError(*err, __LINE__,);
   }

   mvdens_chdump("fishTinv", fish, err);
   forwardError(*err, __LINE__,);

   if (mringstate->maxF==FoM || mringstate->maxF==FoMnocor)  {
      /* Figure of merit = error ellipse area */
      res = figure_of_merit(fish, mringstate->maxF, err);
      forwardError(*err, __LINE__,);
   } else {
      /* Trace */
      for (n=0,res=0.0; n<fish->ndim; n++) res += fish->std[n*fish->ndim+n];
      if (mringstate->maxF==minusTrFm1) res = -res;
      else if (mringstate->maxF==inverseTrFm1) res = 1.0/res;
   }

   mvdens_free(&fish);

   fflush(stdout);
   fprintf(stdout, "%g\n", res);


   /* Print covariance to file */
   if (mringstate->bintype==0) {
      /* Linear bins */
      dtheta = (mringstate->THETA_MAX-mringstate->theta_max_start)/((double)mringstate->Ntheta-1.0);
      dlogtheta = 0.0;
   } else {
      /* Logarithmic bins */
      dtheta = 0.0;
      dlogtheta = (log(mringstate->THETA_MAX)-log(mringstate->theta_max_start))/(mringstate->Ntheta-1.0);
   }

   F = fopen_err("cov_Z", "w", err);
   forwardError(*err, __LINE__,);

   for (i=0,theta_i=mringstate->theta_max_start; i<mringstate->Ntheta; i++) {
      for (j=0,theta_j=mringstate->theta_max_start; j<mringstate->Ntheta; j++) {

	 theta_max_ij[0] = theta_i;
	 theta_max_ij[1] = theta_j;
	 if (eta<=0) {
	    theta_min_ij[0] = theta_min_ij[1] = mringstate->THETA_MIN;
	 } else {
	    theta_min_ij[0] = theta_max_ij[0] * eta;
	    theta_min_ij[1] = theta_max_ij[1] * eta;
	 }

	 //printf("% f % f\n", theta_min_ij[0]/theta_max_ij[0], theta_min_ij[1]/theta_max_ij[1]);

	 cov = cov_RR_Z(theta_min_ij, theta_max_ij, mringstate->cov_xi->theta,
			mringstate->cov_xi->cov, mringstate->cov_xi->n, err);
	 forwardError(*err, __LINE__,);
	 fprintf(F, "%f %f % g\n", theta_max_ij[0]/arcmin, theta_max_ij[1]/arcmin, cov);


	 if (mringstate->bintype==0) theta_j += dtheta;
	 else theta_j *= exp(dlogtheta);

      }

      if (mringstate->bintype==0) theta_i += dtheta;
      else theta_i *= exp(dlogtheta);
   }
   fclose(F);

   free_parameters_lens(&lens);
}

void fisherMap(const char cname[])
{
   error *myerr = NULL, **err;
   int n, diag_only, N;
   double res, theta[2], dtheta, dlogtheta, cov;
   int i, j;
   cosmo_lens *lens;
   config_max config;
   mring_state *mringstate=NULL;
   mvdens *fish;
   FILE *F;


   err = &myerr;

   /* N unused here */
   read_config_and_set_mring(cname, &config, &mringstate, &N, err);
   exitOnError(*err, stderr);

   lens = set_cosmological_parameters_to_default_lens(err);
   exitOnError(*err, stderr);

   diag_only = mringstate->ltype==fisher_var ? 1 : 0;

   fish = fisher_matrix_Map2(lens, mringstate->theta_max_start, mringstate->THETA_MAX,
			     mringstate->Ntheta, mringstate->bintype, mringstate->par, mringstate->npar,
			     diag_only, mringstate->cov_xi, err);
   exitOnError(*err, stderr);

   mvdens_chdump("fishMap", fish, err);
   exitOnError(*err, stderr);

   if (mringstate->maxF==minusTrFm1 || mringstate->maxF==inverseTrFm1 ||
       mringstate->maxF==FoM || mringstate->maxF==FoMnocor) {
      sm2_inverse(fish->std, fish->ndim, err);
      exitOnError(*err, stderr);
   }

   mvdens_chdump("fishMapinv", fish, err);
   exitOnError(*err, stderr);

   if (mringstate->maxF==FoM || mringstate->maxF==FoMnocor)  {
      /* Figure of merit = error ellipse area */
      res = figure_of_merit(fish, mringstate->maxF, err);
      exitOnError(*err, stderr);
   } else {
      /* Trace */
      for (n=0,res=0.0; n<fish->ndim; n++) res += fish->std[n*fish->ndim+n];
      if (mringstate->maxF==minusTrFm1) res = -res;
      else if (mringstate->maxF==inverseTrFm1) res = 1.0/res;
   }

   printf("%g\n", res);

   /* Print covariance to file */
   if (mringstate->bintype==0) {
      /* Linear bins */
      dtheta = (mringstate->THETA_MAX-mringstate->theta_max_start)/((double)mringstate->Ntheta-1.0);
      dlogtheta = 0.0;
   } else {
      /* Logarithmic bins */
      dtheta = 0.0;
      dlogtheta = (log(mringstate->THETA_MAX)-log(mringstate->theta_max_start))/(mringstate->Ntheta-1.0);
   }

   F = fopen_err("cov_Map2", "w", err);
   exitOnError(*err, stderr);      
   for (i=0,theta[0]=mringstate->theta_max_start; i<mringstate->Ntheta; i++) {
      for (j=0,theta[1]=mringstate->theta_max_start; j<mringstate->Ntheta; j++) {

	 cov = cov_Map2(theta, mringstate->cov_xi->theta, mringstate->cov_xi->cov, mringstate->cov_xi->n, err);
	 fprintf(F, "%f %f % g\n", theta[0]/arcmin, theta[1]/arcmin, cov);
	 exitOnError(*err, stderr);      

	 if (mringstate->bintype==0) theta[1] += dtheta;
	 else theta[1] *= exp(dlogtheta);
      }
      if (mringstate->bintype==0) theta[0] += dtheta;
      else theta[0] *= exp(dlogtheta);
   }
   fclose(F);

   mvdens_free(&fish);
   free_parameters_lens(&lens);
}

void fisher_random_a(const char cname[],  const char aname[], error **err)
{
   int i, n, N=0, diag_only;
   double *a, eta, A, B;
   double realA, realB, res, x, tp, tm, y;
   cosmo_lens *lens;
   config_max config;
   mring_state *mringstate=NULL;
   mvdens *fish;
   FILE *F;
   double ran, frac=0.01;
   struct tms time;
   clock_t cl;
   unsigned u;

   read_config_and_set_mring(cname, &config, &mringstate, &N, err);
   forwardError(*err, __LINE__,);

   if (!(mringstate->ltype==fisher_var || mringstate->ltype==fisher_covar)) {
      fprintf(stderr, "Not the correct ltype\n");
      exit(1);
   }

   a = read_a(aname, mringstate->K, N, err);      forwardError(*err, __LINE__,);

   /* K constraints, N total coefficients (including constraints) */
   if (mringstate->THETA_MIN>=0) {
      realA = (mringstate->THETA_MAX - mringstate->THETA_MIN)/2.0;
      realB = (mringstate->THETA_MAX + mringstate->THETA_MIN)/2.0;
   } else {
      realA = (mringstate->THETA_MAX - mringstate->THETA_MAX*mringstate->eta)/2.0;
      realB = (mringstate->THETA_MAX + mringstate->THETA_MAX*mringstate->eta)/2.0;
   }
   if (mringstate->combinef!=all_sc) {
      A = realA;
      B = realB;
      eta = mringstate->eta;
   } else {
      /* General optimisation, independent of scales, set A,B<0 */
      A = B = -1.0;
      eta = -1.0;
   }

   //printf("^^^ %d %d %g %g %d  %f %f\n", N, mringstate->K, A/arcmin, B/arcmin, mringstate->poly,
   //  mringstate->THETA_MIN/arcmin, mringstate->THETA_MAX/arcmin);
   fill_constraint(a, N, mringstate->K, A, B, mringstate->poly, 0, err);
   forwardError(*err, __LINE__,);

   
   /*  for (n=0; n<N; n++) {
       fprintf(stderr, "a[%d] = % .10f; ", n, a[n]);
       }
       fprintf(stderr, "\n");
   */

   /* random init */
   do {
      cl = times(&time);
      u = (unsigned)cl;
      u = (u%10)*(u%100)*(u%1000);
   } while (u==0);
   srand(u);

   for (i=mringstate->K;i<N;i++){
      ran=rand()/(double)RAND_MAX;
      a[i] = a[i] + (2*ran-1)*frac;
    }
 
   fill_constraint(a, N, mringstate->K, A, B, mringstate->poly, 0, err);
   exitOnError(*err, stderr);
   normalise(a, N, mringstate->poly, err);
   exitOnError(*err, stderr);

   lens = set_cosmological_parameters_to_default_lens(err);
   forwardError(*err, __LINE__,);

   //printf("THETA_MIN eta A B theta_max_start = %g' %g  %g' %g'  %g'\n",
   //	  mringstate->THETA_MIN/arcmin, eta, A/arcmin, B/arcmin, mringstate->theta_max_start/arcmin);

   diag_only = mringstate->ltype==fisher_var ? 1 : 0;
   fish = fisher_matrix_mring(lens, mringstate->THETA_MIN, mringstate->theta_max_start, mringstate->THETA_MAX,
			      mringstate->Ntheta, mringstate->bintype, a, N, mringstate->poly, mringstate->par,
			      mringstate->npar, diag_only, mringstate->cov_xi, eta, NULL, err);
   forwardError(*err, __LINE__,);

   mvdens_chdump("fishT", fish, err);
   forwardError(*err, __LINE__,);
   
   if (mringstate->maxF==minusTrFm1 || mringstate->maxF==inverseTrFm1 ||
       mringstate->maxF==FoM || mringstate->maxF==FoMnocor) {
      sm2_inverse(fish->std, fish->ndim, err);
      forwardError(*err, __LINE__,);
   }

   mvdens_chdump("fishTinv", fish, err);
   forwardError(*err, __LINE__,);

   if (mringstate->maxF==FoM || mringstate->maxF==FoMnocor)  {
      /* Figure of merit = error ellipse area */
      res = figure_of_merit(fish, mringstate->maxF, err);
      forwardError(*err, __LINE__,);
   } else {
      /* Trace */
      for (n=0,res=0.0; n<fish->ndim; n++) res += fish->std[n*fish->ndim+n];
      if (mringstate->maxF==minusTrFm1) res = -res;
      else if (mringstate->maxF==inverseTrFm1) res = 1.0/res;
   }

   fflush(stdout);
   F = fopen("fisher", "w");
   fprintf(F, "%g\n", res);
   fclose(F);



   F = fopen("T+", "w");
   for (x=-1.0; x<=1.001; x+=0.01) {
      tp = Tp(x, a, N, mringstate->poly, err);
      forwardError(*err, __LINE__,);
      tm = Tm(x, a, N, mringstate->poly, realB/realA, err);
      forwardError(*err, __LINE__,);
      y = eta + (1.0-eta)/2.0*(1.0+x);
      fprintf(F, "%f % f % g % g\n", y, x, tp, tm);
   }
   fclose(F);

   mvdens_free(&fish);
   free_parameters_lens(&lens);
   free(a);
}


void plot_cov_RR(const char cname[],  const char aname[], error **err)
{
   int n, N=0, i, j;
   double *a, eta, A, B;
   double theta_min_ij[2], theta_max_ij[2], theta_i, theta_j, cov, dtheta, dlogtheta;
   double realA, realB, nice_plot_fac;
   config_max config;
   mring_state *mringstate=NULL;
   FILE *F;


   read_config_and_set_mring(cname, &config, &mringstate, &N, err);
   forwardError(*err, __LINE__,);

   if (!(mringstate->ltype==fisher_var || mringstate->ltype==fisher_covar)) {
      fprintf(stderr, "Not the correct ltype\n");
      exit(1);
   }

   a = read_a(aname, mringstate->K, N, err);      forwardError(*err, __LINE__,);

   /* K constraints, N total coefficients (including constraints) */
   if (mringstate->THETA_MIN>=0) {
      realA = (mringstate->THETA_MAX - mringstate->THETA_MIN)/2.0;
      realB = (mringstate->THETA_MAX + mringstate->THETA_MIN)/2.0;
   } else {
      realA = (mringstate->THETA_MAX - mringstate->THETA_MAX*mringstate->eta)/2.0;
      realB = (mringstate->THETA_MAX + mringstate->THETA_MAX*mringstate->eta)/2.0;
   }
   if (mringstate->combinef!=all_sc) {
      A = realA;
      B = realB;
      eta = mringstate->eta;
   } else {
      /* General optimisation, independent of scales, set A,B<0 */
      A = B = -1.0;
      eta = -1.0;
   }

   //printf("^^^ %d %d %g %g %d  %f %f\n", N, mringstate->K, A/arcmin, B/arcmin, mringstate->poly,
   //  mringstate->THETA_MIN/arcmin, mringstate->THETA_MAX/arcmin);
   fill_constraint(a, N, mringstate->K, A, B, mringstate->poly, 0, err);
   forwardError(*err, __LINE__,);

   //normalise(a, N, mringstate->poly, err);
   //forwardError(*err, __LINE__,);


   for (n=0; n<N; n++) {
      fprintf(stderr, "a[%d] = % .10f; ", n, a[n]);
   }
   fprintf(stderr, "\n");

   /* Print covariance to file */
   if (mringstate->bintype==0) {
      /* Linear bins */
      dtheta = (mringstate->THETA_MAX-mringstate->theta_max_start)/((double)mringstate->Ntheta-1.0);
      dlogtheta = 0.0;
   } else {
      /* Logarithmic bins */
      dtheta = 0.0;
      dlogtheta = (log(mringstate->THETA_MAX)-log(mringstate->theta_max_start))/(mringstate->Ntheta-1.0);
   }

   F = fopen_err("cov_RR", "w", err);
   forwardError(*err, __LINE__,);

   nice_plot_fac = 1.02;

   for (i=0,theta_i=mringstate->theta_max_start; i<mringstate->Ntheta; i++) {
      for (j=0,theta_j=mringstate->theta_max_start; j<mringstate->Ntheta; j++) {

	 theta_max_ij[0] = theta_i;
	 theta_max_ij[1] = theta_j;
	 if (eta<=0) {
	    theta_min_ij[0] = theta_min_ij[1] = mringstate->THETA_MIN;
	 } else {
	    theta_min_ij[0] = theta_max_ij[0] * eta;
	    theta_min_ij[1] = theta_max_ij[1] * eta;
	 }

	 //printf("% f % f\n", theta_min_ij[0]/theta_max_ij[0], theta_min_ij[1]/theta_max_ij[1]);

	 cov = cov_RR(theta_min_ij, theta_max_ij, a, N, mringstate->poly, mringstate->cov_xi->theta,
		      mringstate->cov_xi->cov, mringstate->cov_xi->n,
		      outer, -1, -1, nice_plot_fac, err);
	 forwardError(*err, __LINE__,);
	 fprintf(F, "%f %f % g\n", theta_max_ij[0]/arcmin, theta_max_ij[1]/arcmin, cov);


	 if (mringstate->bintype==0) theta_j += dtheta;
	 else theta_j *= exp(dlogtheta);

      }

      if (mringstate->bintype==0) theta_i += dtheta;
      else theta_i *= exp(dlogtheta);
   }
   fclose(F);
}


void usage()
{
   fprintf(stderr, "Usage: mringtest [task] [config] [list_a] [...]\n");
   fprintf(stderr, "(Defaults: config_test, list_a)\n");
   fprintf(stderr, "   0 S_N_T           config list_a      ");
   fprintf(stderr, "     read coefficients from list_a, output S/N\n");
   fprintf(stderr, "   1 S_N_Z_comb_none config             ");
   fprintf(stderr, "     S/N for Z+\n");
   fprintf(stderr, "   3 S_N_sum         config             ");
   fprintf(stderr, "     Sum of S/N for Z+\n");
   fprintf(stderr, "   4 S_N_Map         config             ");
   fprintf(stderr, "     S/N for Map^2\n");
   fprintf(stderr, "   5 fisherT         config list_a      ");
   fprintf(stderr, "     fisher for T+\n");
   fprintf(stderr, "   6 extend_Tp_zero  config list_a psi_1");
   fprintf(stderr, "     Extend T+ to larger psi, fill with zero\n");
   fprintf(stderr, "   7 fisherZ         config             ");
   fprintf(stderr, "     fisher for Z+\n");
   fprintf(stderr, "   8 fisherMap       config             ");
   fprintf(stderr, "     fisher for Map^2\n");
   fprintf(stderr, "   9 check_EB        config list_a      ");
   fprintf(stderr, "     E-/B-moe decomposition using T+, T-\n");
   fprintf(stderr, "   10 S_N_random_a       config list_a      ");
   fprintf(stderr, "     SN for randomisation of an \n");
   fprintf(stderr, "   11 fisher_random_a       config list_a      ");
   fprintf(stderr, "      fisher for randomisation of an \n");
   fprintf(stderr, "   11 plot_cov_RR    config list_a      ");
   fprintf(stderr, "      Covariance for RE\n");
}

int main(int argc, char *argv[])
{
   int task;
   char cname[128], aname[128];
   double psi1 = 0.0;
   error *myerr = NULL, **err;

   err = &myerr;

   /* Default options */
   strcpy(cname, "config_test");
   strcpy(aname, "list_a");


   switch (argc) {
      case 5 :
	 psi1 = atof(argv[4]);
      case 4 :
	 strcpy(aname, argv[3]);
      case 3 :
	 strcpy(cname, argv[2]);
      case 2 : 
	 task = atoi(argv[1]);
	 break;
      default :
	 usage();
	 exit(1);
   }

   fprintf(stderr, "task=%d configname=%s aname=%s\n", task, cname, aname);


   random_init("init");
   //    srand(31416); fprintf(stderr, "*** NO RANDOM INIT ***\n");


   switch (task) {
      case 0 :
	 /* Read list_a, calculate S/N *
	  * Not for all_sc_sum!!!      *
	  * Old name compare_T_R       */
	 S_N_T(cname, aname);
	 break;
      case 1 :
	 /* S/N for Z+ */
	 S_N_Z_comb_none(cname);
	 break;
      case 3 :
	 S_N_sum(cname);
	 break;
      case 4 :
	 S_N_Map(cname);
	 break;
      case 5 :
	 fisherT(cname, aname, err);
	 exitOnError(*err, stderr);
	 break;
      case 6 :
	 extend_Tp_zero(cname, aname, psi1*arcmin);
	 break;
      case 7 :
	 fisherZ(cname, err);
	 exitOnError(*err, stderr);
	 break;
      case 8 :
	 fisherMap(cname);
	 break;
      case 9 :
	 check_EB(cname, aname);
	 break;
      case 10 :
	 S_N_random_a(cname, aname);
	 break;
      case 11 :
	 fisher_random_a(cname, aname, err);
	 exitOnError(*err, stderr);
      case 12 :
	 plot_cov_RR(cname, aname, err);
	 exitOnError(*err, stderr);
	 break;

      default :
	 usage();
	 exit(2);
   }
  
   return 0;
}

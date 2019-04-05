#include <math.h>

#include "param.h"
#include "xcorr.h"


#define NR_END 1
#define FREE_ARG char*

#ifndef DSQR
static double darg __attribute__((unused));
#define DSQR(a) ((darg=(a))==0.0 ? 0.0 : darg*darg)
#endif

double *vector(long nl, long nh, error **err)
/* allocate a double vector with subscript range v[nl..nh] */
{
   double *v;

   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   testErrorRet(!v, xc_allocate, "allocation failure in vector()", *err, __LINE__, NULL);
   return v-nl+NR_END;
}

void free_vector(double *v)
/* free a double vector allocated with vector() */
{
   free((FREE_ARG) (v-NR_END));
}

double **matrix(long nrl, long nrh, long ncl, long nch, error **err)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double **m;

   /* allocate pointers to rows */
   m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
   testErrorRet(!m, xc_allocate, "allocation failure 1 in matrix()", *err, __LINE__, NULL);
   m += NR_END;
   m -= nrl;

   /* allocate rows and set pointers to them */
   m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
   testErrorRet(!m[nrl], xc_allocate, "allocation failure 2 in matrix()", *err, __LINE__, NULL);
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   /* return pointer to array of pointers to rows */
   return m;
}

void free_matrix(double **m)
{
   free((FREE_ARG) (m[0]-NR_END));
   free((FREE_ARG) (m-NR_END));
}
#undef FREE_ARG
#undef NR_END

void matrix_print(const double **a, int n, FILE *F)
{
   int i, j;

   for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
	 fprintf(F, "% .4f ", a[i][j]);
      }
      fprintf(F, "\n");
   }
}

/* Some of those functions exist in maths.c for matrices as single pointers */

double **matrix_copy(const double **M, int N, error **err)
{
   int i, j;
   double **C;

   C = matrix(0, N-1, 0, N-1, err);
   forwardError(*err, __LINE__, NULL);
   for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
	 C[i][j] = M[i][j];
      }
   }
   return C;
}

/* ============================================================ *
 * (See NR p.46 - 49.)						*
 * Inverts the NxN matrix C and stores the inverse in C         *
 * (old matrix gets destroyed). Returns det C.			*
 * The index goes from 0 to N-1!				*
 * ============================================================ */

double matrix_inverse(double **C, int N, error **err)
{
   int i, j;
   double *col, det, **Cinv;
   int *p;

   Cinv = (double**)malloc(N*sizeof(double));
   testErrorRet(Cinv==NULL, xc_allocate, "Cannot allocate memory", *err, __LINE__, 0);
   for (i=0; i<N; i++) {
      Cinv[i]   = vector(0, N-1, err);   forwardError(*err, __LINE__, 0);
   }
   col = vector(0, N-1, err);            forwardError(*err, __LINE__, 0);
   p = (int*)malloc(N*sizeof(int));
   testErrorRet(p==NULL, xc_allocate, "Cannot allocate memory", *err, __LINE__, 0);

   matrix_ludcmp(C, N, p, &det, err);    forwardError(*err, __LINE__, 0);
   for (j=0; j<N; j++) {
      det *= C[j][j];
   }

   for (j=0; j<N; j++) {
      for (i=0; i<N; i++) {
	 col[i] = 0.0;
      }
      col[j] = 1.0;
      matrix_lubksb(C, N, p, col, err);  forwardError(*err, __LINE__, 0);
      for (i=0; i<N; i++) {
	 Cinv[i][j] = col[i];
      }
   }
   for (j=0; j<N; j++) {
      for (i=0; i<N; i++) {
	 C[i][j] = Cinv[i][j];
      }
   }
   for (i=0; i<N; i++) {
      free_vector(Cinv[i]);
   }
   free((void*)Cinv);

   return det;
}

#define EEPS 1.0e-15
void matrix_lubksb(double **a, int n, int *indx, double b[], error **err)
{
   int i,ii=-1,ip,j;
   double sum;

   for (i=0;i<n;i++) {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii!=-1)
	for (j=ii;j<=i-1;j++) {
	   sum -= a[i][j]*b[j];
	}
      else if (sum) ii=i;
      b[i]=sum;
   }
   for (i=n-1;i>=0;i--) {
      sum=b[i];
      for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
      testErrorRet(fabs(a[i][i])<EEPS, xc_infnan, "Division by zero", *err, __LINE__,);
      b[i]=sum/a[i][i];
   }
}
#undef EEPS

#define TINY 1.0e-20
void matrix_ludcmp(double **a, int n, int *indx, double *d, error **err)
{
   int i,imax=-1,j,k;
   double big,dum,sum,temp;
   double *vv;

   vv=vector(0,n-1, err);        forwardError(*err, __LINE__,);
   *d=1.0;
   for (i=0;i<n;i++) {
      big=0.0;
      for (j=0;j<n;j++)
	if ((temp=fabs(a[i][j])) > big) big=temp;
      testErrorRet(big == 0.0, xc_singular, "Singular matrix in routine ludcmp", *err, __LINE__,);
      vv[i]=1.0/big;
   }
   for (j=0;j<n;j++) {
      for (i=0;i<j;i++) {
	 sum=a[i][j];
	 for (k=0;k<i;k++)
	   sum -= a[i][k]*a[k][j];
	 a[i][j]=sum;
      }
      big=0.0;
      for (i=j;i<n;i++) {
	 sum=a[i][j];
	 for (k=0;k<j;k++)
	   sum -= a[i][k]*a[k][j];
	 a[i][j]=sum;
	 if ( (dum=vv[i]*fabs(sum)) >= big) {
	    big=dum;
	    imax=i;
	 }
      }
      if (j != imax) {
	 for (k=0;k<n;k++) {
	    dum=a[imax][k];
	    a[imax][k]=a[j][k];
	    a[j][k]=dum;
	 }
	 *d = -(*d);
	 vv[imax]=vv[j];
      }
      indx[j]=imax;
      if (a[j][j] == 0.0) a[j][j]=TINY;
      if (j != n-1) {
	 dum=1.0/(a[j][j]);
	 for (i=j+1;i<n;i++) a[i][j] *= dum;
      }
   }
   free_vector(vv);
}
#undef TINY

double **matrix_multiply(const double **a, const double **b, int N, error **err)
{
   double **c, sum;
   int i, j, k;

   c = matrix(0, N-1, 0, N-1, err);
   forwardError(*err, __LINE__, 0);
   for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
	 for (k=0,sum=0.0; k<N; k++) {
	    sum += a[i][k]*b[k][j];
	 }
	 c[i][j] = sum;
      }
   }

   return c;
}

#define maxdev 1.0e-5
void matrix_unity_test(const double **a, const double **ainv, int N, error **err)
{
   int i, j;
   double max1, max2;

   double **unit;

   unit = matrix_multiply(a, ainv, N, err);
   forwardError(*err, __LINE__,);

   for (i=0,max1=-1.0,max2=-1.0; i<N; i++) {
      if (fabs(unit[i][i]-1.0)>max1) max1 = fabs(unit[i][i]-1.0);
      for (j=0; j<N; j++) {
	 if (i!=j && fabs(unit[i][j])>max2) max2 = fabs(unit[i][j]);
      }
   }

   if (max1>maxdev || max2>maxdev) {
      fprintf(stderr, "max |diag-1| = %.2e, max |offdiag| = %.2e\n", max1, max2);
      *err = addError(xc_unity, "unit matrix test failed", *err, __LINE__);
   }
   free_matrix(unit);
}

functions_wrapper_t *init_functions_Xcorr(error **err)
{
   functions_wrapper_t *init_func;

   fprintf(stderr, "init_functions_Xcorr\n");

   init_func = init_functions_wrapper(read_from_config_Xcorr, init_Xcorr, likeli_Xcorr, NULL, NULL, print_Xcorr, err);
   forwardError(*err, __LINE__, NULL);

   return init_func;
}
 
void read_from_config_Xcorr(void **state, FILE *F, error **err)
{
   xcorr_state *xcorrstate;
   config_element c = {0, 0.0, ""};

   xcorrstate = malloc_err(sizeof(xcorr_state), err);  forwardError(*err, __LINE__,);
   CONFIG_READ_S(xcorrstate, datname, s, F, c, err);
   *state = xcorrstate;
}

void init_Xcorr(common_like *like, error **err)
{
   xcorr_state *xcorr;

   xcorr = (xcorr_state*)(like->state);

   testErrorRet(xcorr==NULL, xc_allocate, "Cannot allocate", *err, __LINE__,);

   read_xcorr(xcorr->datname, xcorr, err);
   forwardError(*err, __LINE__,);
}

void read_xcorr(const char *name, xcorr_state *xc, error **err)
{
   int i, j, nread;
   FILE *F;

   F = fopen(name, "r");
   testErrorRet(F==NULL, xc_file, "Could not open file", *err, __LINE__,);
   nread = fscanf(F, "# Nth = %d Mz = %d\n", &xc->Nth, &xc->Mz);
   testErrorRet(nread!=2, xc_io, "data file has not the right format (line 1)", *err, __LINE__,);

   xc->theta   = malloc(sizeof(double)*xc->Nth);
   xc->nbar    = malloc(sizeof(double)*xc->Mz);
   xc->signbar = malloc(sizeof(double)*xc->Mz);
   xc->xi      = matrix(0, xc->Nth-1, 0, xc->Mz*(xc->Mz+1)/2-1, err);
   forwardError(*err, __LINE__,);
   xc->sigxi   = matrix(0, xc->Nth-1, 0, xc->Mz*(xc->Mz+1)/2-1, err);
   forwardError(*err, __LINE__,);

   nread = fscanf(F, "# ");
   for (j=0; j<xc->Mz; j++) {
      nread = fscanf(F, "%lg %lg ", xc->nbar+j, xc->signbar+j);
      testErrorRet(nread!=2, xc_io, "data file has not the right format (line 2)", *err, __LINE__,);
   }

   for (i=0; i<xc->Nth; i++) {
      nread = fscanf(F, "%lg ", xc->theta+i);
      testErrorRet(nread!=1, xc_io, "data file has not the right format (line 3+, column 1)", *err, __LINE__,);
      for (j=0; j<xc->Mz*(xc->Mz+1)/2; j++) {
	 nread = fscanf(F, "%lg %lg", &(xc->xi[i][j]), &(xc->sigxi[i][j]));
	 //fprintf(stderr, "%d %d %d    % g % g\n", j, i, nread, (xc->xi[i][j]), (xc->sigxi[i][j]));
         testErrorRet(nread!=2, xc_io, "data file has not the right format (line 3+, column 2+)", *err, __LINE__,);
      }
   }

}

xcorr_model *init_xcorr_model(int Nth, int Mz, error **err)
{
   xcorr_model *model;

   model         = malloc(sizeof(xcorr_model));            forwardError(*err, __LINE__, NULL);

   model->eta    = matrix(0, Mz-1, 0, Mz-1, err);          forwardError(*err, __LINE__, NULL);

   model->A      = malloc(sizeof(double)*Nth);
   testErrorRet(model->A==NULL, xc_allocate, "Cannot allocate", *err, __LINE__, NULL);

   model->gamma  = malloc(sizeof(double)*Nth);
   testErrorRet(model->gamma==NULL, xc_allocate, "Cannot allocate", *err, __LINE__, NULL);

   model->Nth = Nth;
   model->Mz  = Mz;

   return model;
}

void free_xcorr_model(xcorr_model *model)
{
   free_matrix(model->eta);
   free_matrix(model->etainv);
   free(model->A);
   free(model->gamma);
}

xcorr_model *params_to_xcorr_model(const xcorr_state *xcorr, const double *params, error **err)
{
   xcorr_model *model;
   int Nth, Mz;
   /*int i, j, k;
     double sum;*/

   Nth = xcorr->Nth;
   Mz  = xcorr->Mz;

   model = init_xcorr_model(Nth, Mz, err);
   forwardError(*err, __LINE__, NULL);

   /*
   switch (xcorr->param) {

      case param_pl :

	 for (j=0,k=0; j<Mz; j++) {
	    for (i=0,sum=0; i<Mz; i++) {
	       if (i==j) continue;
	       model->eta[i][j] = params[k++];
	       sum += model->eta[i][j];
	    }
	    model->eta[j][j] = 1.0-sum;
	 }
	 for (i=0; i<Mz; i++) {
	    model->A[i] = params[k++];
	 }
	 for (i=0; i<Mz; i++) {
	    model->gamma[i] = params[k++];
	 }
	 break;

      case param_pl_gfix :

	 for (j=0,k=0; j<Mz; j++) {
	    for (i=0,sum=0.0; i<Mz; i++) {
	       if (i==j) continue;
	       model->eta[i][j] = params[k++];
	       sum += model->eta[i][j];
	    }
	    model->eta[j][j] = 1.0-sum;
	 }
	 for (i=0; i<Mz; i++) {
	    model->A[i] = params[k++];
	 }
	 for (i=0; i<Mz; i++) {
	    model->gamma[i] = 0.18;
	 }
	 break;

      default :

	 *err = addError(xc_unknown, "unknown parametrization", *err, __LINE__);
	 return NULL;

   }
   */

   model->etainv = matrix_copy((const double**)model->eta, Mz, err);
   forwardError(*err, __LINE__, NULL);
   matrix_inverse(model->etainv, Mz, err);                forwardError(*err, __LINE__, NULL);
   matrix_unity_test((const double**)model->eta, (const double**)model->etainv, Mz, err);
   forwardError(*err, __LINE__, NULL);

   return model;
}

double likeli_Xcorr(common_like *like, const double *params, error **err)
{
   xcorr_state *xcorr;
   xcorr_model *model;
   double res;

   xcorr = (xcorr_state*)(like->state);

   model = params_to_xcorr_model(xcorr, params, err);   forwardError(*err, __LINE__, 0);
   res = xcorr_chi2(xcorr, model, err);                 forwardError(*err, __LINE__, 0);

   if (model) free_xcorr_model(model);

   return res;
}

double w_true(double theta, double A, double gamma, error **err)
{
   testErrorRet(theta<0, xc_infnan, "Division by zero (or worse)", *err, __LINE__, 0);
   return A*pow(theta, -gamma);
}

double w_photoz(double theta, int i, int j, const xcorr_model *model, const double *nbar, error **err)
{
   int k, l, Mz;
   double sum, x, y, wtrue;

   Mz = model->Mz;
   testErrorRet(i<0 || j<0 || i>=Mz || j>=Mz, xc_range, "Redshift bin indices out of range",
		*err, __LINE__, 0.0);

   for (k=0,sum=0.0; k<Mz; k++) {
      wtrue = w_true(theta, model->A[k], model->gamma[k], err);
      forwardError(*err, __LINE__, 0);

      x = model->eta[i][k]*model->eta[j][k];
      for (l=0,y=0.0; l<Mz; l++) {
	 y += model->etainv[k][l]*nbar[l];
      }
      sum += x*y*y*wtrue;
   }

   x = nbar[i]*nbar[j];
   testErrorRet(x<=0, xc_infnan, "Division by zero (or worse)", *err, __LINE__, 0);
   return sum/x;
}

double xcorr_chi2(const xcorr_state *xcorr, const xcorr_model *model, error **err)
{
   int i, j, k, m, Mz;
   double chi2, wph;

   Mz = model->Mz;
   testErrorRet(Mz!=xcorr->Mz, xc_inconsistent, "Number of z-bins in model and state different",
		*err, __LINE__, 0);

   //matrix_print(model->eta, model->Mz, stdout);
   //matrix_print(model->etainv, model->Mz, stdout);
   for (i=0,m=0,chi2=0.0; i<Mz; i++) {
      for (j=0; j<=i; j++,m++) {
	 for (k=0; k<xcorr->Nth; k++) {
	    wph = w_photoz(xcorr->theta[k], i, j, model, xcorr->nbar, err);
	    forwardError(*err, __LINE__, 0);

	    /* TODO: error on nbar */
	    chi2 += DSQR((wph - xcorr->xi[k][m])/xcorr->sigxi[k][m]);
	    //if (i==Mz-1 && k==0)  printf("chi2: model - obs, err =  % .3f % .3f % .3f   %g\n", wph,
	    //		      xcorr->xi[k][m], xcorr->sigxi[k][m], chi2);
	 }
      }
   }

   chi2 = -0.5*chi2;
   return chi2;
}

void print_Xcorr(FILE *F, void *state, error **err)
{
   xcorr_state *xc;
   FILE *rhere;
   int i, k, Mz;

   if (F==NULL) rhere = stdin;
   else rhere = F;

   xc = (xcorr_state*)state;

   fprintf(rhere, "Nth = %d, Mz = %d\n", xc->Nth, xc->Mz);
   Mz = xc->Mz;
   fprintf(F, "nbar, signbar:\n");
   for (i=0; i<Mz; i++) {
      fprintf(rhere, "%g %g   ", xc->nbar[i], xc->signbar[i]);
   }
   fprintf(rhere, "\ntheta xi, sigxi:\n");
   for (k=0; k<xc->Nth; k++) {
      fprintf(rhere, "%.3f   ", xc->theta[k]);
      for (i=0; i<Mz*(Mz+1)/2; i++) {
	 fprintf(rhere, "%g %g   ", xc->xi[k][i], xc->sigxi[k][i]);
      }
      fprintf(rhere, "\n");
   }
}

void out_xcorr_model(const xcorr_model *model, const double *theta, FILE *F, error **err)
{
   int i, k;
   FILE *out;

   if (F==NULL) out = stdout;
   else out = F;

   fprintf(out, "model:\n");
   fprintf(out, "======\n");
   fprintf(out, "eta:\n");
   matrix_print((const double**)model->eta, model->Mz, out);
   fprintf(out, "etainv:\n");
   matrix_print((const double **)model->etainv, model->Mz, out);
   fprintf(out, "A:\n");
   for (i=0; i<model->Mz; i++) {
      fprintf(out, "% 6.2f ", model->A[i]);
   }
   fprintf(out, "\ngamma:\n");
   for (i=0; i<model->Mz; i++) {
      fprintf(out, "% .3f ", model->gamma[i]);
   }

   fprintf(out, "\n# theta w_true(theta)\n");
   for (k=0; k<model->Nth; k++) {
      fprintf(out, "%.3f ", theta[k]);
      for (i=0; i<model->Mz; i++) {
	 fprintf(out, "% g ", w_true(theta[k], model->A[i], model->gamma[i], err));
	 forwardError(*err, __LINE__,);
      }
      fprintf(out, "\n");
   }
   fprintf(out, "======\n");
}

/* ============================================================ *
 * mkmax.c							*
 * Martin Kilbinger 2008, 2009					*
 * ============================================================ */

#include "mkmax.h"

/* Posterior log pdf for the cg maximum search */
double post_cg(config_base *config, const double *x, const parabox *pb, int quiet, error **err)
{
   double val;
   int inside;
  
   inside = isinBox(pb, x, err);      forwardError(*err, __LINE__, 0);
   if (inside==0) {
      val = MC_LOGBIG;
   } else {
     /*val = -log P = -[-0.5*chi2 + log(prior)] */
     val = -posterior_log_pdf_common(config, x, err);
   }

   if (getErrorValue(*err)!=noErr) {
      val = MC_LOGBIG;

      if (getErrorValue(*err)==tls_cosmo_par) {
	 forwardError(*err, __LINE__, val);
      } else {
	 do { 	   /* Dummy loop */
	    ParameterErrorVerb(*err, x, quiet, config->npar);
	 } while (0);
      }
   }

   if (! quiet) {
      fprintf(stderr, "post_cg: log P = % .6f   par = ", -val);
      print_parameter(stderr, config->npar, x);
   }
  
   return val;
}

/* Numerical derivatie of the posterior */
#define eps 1.0e-20
void dposterior_log_pdf(config_base *config, const double *x, double *df, error **err)
{
   int i;
   double h, *param, a, b;

   param     = sm2_vector(0, config->npar-1, err); forwardError(*err, __LINE__,);

   for (i=0; i<config->npar; i++) {
      param[i] = x[i];
   }

   for (i=0; i<config->npar; i++) {
      h = fh*(config->max[i] - config->min[i]);                 /* fh times box size */
      testErrorRet(h<eps, ce_infnan, "h too small", *err, __LINE__,);

      param[i] += h;
      a = posterior_log_pdf_common(config, param, err);
      if (getErrorValue(*err)!=noErr) {
	 purgeError(err);
	 a = MC_LOGBIG;
      }

      param[i] -= 2.0*h;
      b = posterior_log_pdf_common(config, param, err);
      if (getErrorValue(*err)!=noErr) {
	 purgeError(err);
	 b = MC_LOGBIG;
      }

      df[i] = (a-b)/(2.0*h);
      param[i] += h;
   }

   sm2_free_vector(param, 0, config->npar-1);
}
#undef eps

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)
#define FMAX(a,b) (a>b ? a : b)
void mnbrak(config_base *config, double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, 
	    const double *p, const double *xi, const parabox *pb, int quiet, error **err)
{
   double ulim, u, r, q, fu, dum, *xt, tmp;
   int n, i;

   n = config->npar;
   xt = sm2_vector(0, n-1, err);                       forwardError(*err, __LINE__,);
   testErrorRet(!finite(*ax), math_infnan, "Parameter *ax is not finite", *err, __LINE__,);
   for (i=0; i<n; i++) {
      xt[i] = p[i] + (*ax)*xi[i];
      testErrorRetVA(!finite(xi[i]), math_infnan, "Parameter xi[%d] is not finite", *err, __LINE__,, i);
      testErrorRetVA(!finite(xt[i]), math_infnan, "Parameter xi[%d] is not finite", *err, __LINE__,, i);
      testErrorRetVA(!finite(p[i]), math_infnan, "Parameter p[%d] is not finite", *err, __LINE__,, i);
   }

   *fa = post_cg(config, xt, pb, quiet, err);     forwardError(*err, __LINE__,);

   for (i=0; i<n; i++) xt[i] = p[i] + *bx*xi[i];
   *fb = post_cg(config, xt, pb, quiet, err);     forwardError(*err, __LINE__,);

   if (*fb > *fa) {
      SHFT(dum,*ax,*bx,dum);
      SHFT(dum,*fb,*fa,dum);
   }
   *cx=(*bx)+GOLD*(*bx-*ax);
   for (i=0; i<n; i++) xt[i] = p[i] + *cx*xi[i];
   *fc = post_cg(config, xt, pb, quiet, err);     forwardError(*err, __LINE__,);
   while (*fb > *fc) {
      r=(*bx-*ax)*(*fb-*fc);
      q=(*bx-*cx)*(*fb-*fa);
      u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
	(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
      ulim=(*bx)+GLIMIT*(*cx-*bx);
      if ((*bx-u)*(u-*cx) > 0.0) {
	 for (i=0; i<n; i++) xt[i] = p[i] + u*xi[i];
	 fu = post_cg(config, xt, pb, quiet, err);   forwardError(*err, __LINE__,);
	 if (fu < *fc) {
	    *ax=(*bx);
	    *bx=u;
	    *fa=(*fb);
	    *fb=fu;
	    goto end;
	 } else if (fu > *fb) {
	    *cx=u;
	    *fc=fu;
	    goto end;
	 }
	 u=(*cx)+GOLD*(*cx-*bx);
	 for (i=0; i<n; i++) xt[i] = p[i] + u*xi[i];
	 fu = post_cg(config, xt, pb, quiet, err);
      } else if ((*cx-u)*(u-ulim) > 0.0) {
	 for (i=0; i<n; i++) xt[i] = p[i] + u*xi[i];
	 fu = post_cg(config, xt, pb, quiet, err);
	 if (fu < *fc) {
	    SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
	    for (i=0; i<n; i++) xt[i] = p[i] + u*xi[i];
	    tmp = post_cg(config, xt, pb, quiet, err);
	    SHFT(*fb,*fc,fu,tmp);
	 }
      } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
	 u=ulim;
	 for (i=0; i<n; i++) xt[i] = p[i] + u*xi[i];
	 fu = post_cg(config, xt, pb, quiet, err);
      } else {
	 u=(*cx)+GOLD*(*cx-*bx);
	 for (i=0; i<n; i++) xt[i] = p[i] + u*xi[i];
	 fu = post_cg(config, xt, pb, quiet, err);
      }
      SHFT(*ax,*bx,*cx,u);
      SHFT(*fa,*fb,*fc,fu);
   }

 end:
   sm2_free_vector(xt, 0, n-1);
   return;
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)
double brent(config_base *config, double ax, double bx, double cx, double tol, double *xmin,
	     const double *pp, const double *xi, const parabox *pb, int quiet, error **err)
{
   int iter, i, n;
   double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm, *xt;
   double e=0.0;

   n  = config->npar;
   xt = sm2_vector(0, n-1, err);			  forwardError(*err, __LINE__, 0);

   a=(ax < cx ? ax : cx);
   b=(ax > cx ? ax : cx);
   x=w=v=bx;

   for (i=0; i<n; i++) xt[i] = pp[i] + x*xi[i];
   fw=fv=fx=post_cg(config, xt, pb, quiet, err);   forwardError(*err, __LINE__, 0);
   for (iter=1;iter<=ITMAX;iter++) {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
      if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
	 *xmin=x;
	 sm2_free_vector(xt, 0, n-1);
	 return fx;
      }
      if (fabs(e) > tol1) {
	 r=(x-w)*(fx-fv);
	 q=(x-v)*(fx-fw);
	 p=(x-v)*q-(x-w)*r;
	 q=2.0*(q-r);
	 if (q > 0.0) p = -p;
	 q=fabs(q);
	 etemp=e;
	 e=d;
	 if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	   d=CGOLD*(e=(x >= xm ? a-x : b-x));
	 else {
	    d=p/q;
	    u=x+d;
	    if (u-a < tol2 || b-u < tol2)
	      d=SIGN(tol1,xm-x);
	 }
      } else {
	 d=CGOLD*(e=(x >= xm ? a-x : b-x));
      }
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      for (i=0; i<n; i++) xt[i] = pp[i] + u*xi[i];
      fu = post_cg(config, xt, pb, quiet, err);      forwardError(*err, __LINE__, 0);
      if (fu <= fx) {
	 if (u >= x) a=x; else b=x;
	 SHFT(v,w,x,u);
	 SHFT(fv,fw,fx,fu);
	   } else {
	      if (u < x) a=u; else b=u;
	      if (fu <= fw || w == x) {
		 v=w;
		 w=u;
		 fv=fw;
		 fw=fu;
	      } else if (fu <= fv || v == x || v == w) {
		 v=u;
		 fv=fu;
	      }
	   }
   }

   sm2_free_vector(xt, 0, n-1);
   *err = addError(mcmc_tooManySteps, "too many steps", *err, __LINE__);
   return 0.0;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT


#define TOL 2.0e-4
void linmin(config_base *config, double p[], double xi[], double *fret,
	    const parabox *pb, int quiet, error **err)
{
   int j, n;
   double xx,xmin,fx,fb,fa,bx,ax;

   n  = config->npar;
   ax = 0.0;
   xx = 1.0;
   mnbrak(config,&ax,&xx,&bx,&fa,&fx,&fb,p,xi,pb,quiet,err);   forwardError(*err, __LINE__,);
   *fret=brent(config,ax,xx,bx,TOL,&xmin,p,xi,pb,quiet,err);   forwardError(*err, __LINE__,);
   for (j=0;j<n;j++) {
      xi[j] *= xmin;
      p[j]  += xi[j];
   }
}
#undef TOL


/* ============================================================ *
 * NR minimization routine, conjugate-gradient method		*
 * ============================================================ */
#define ITMAX 200
#define ZEPS 1.0e-10
void frprmn(config_base *config, double p[], double ftol, int *iter, double *fret,
	    const parabox *pb, int quiet, error **err)
{
   int j, its, n;
   double gg,gam,fp,dgg;
   double *g,*h,*xi;

   n = config->npar;
   g = sm2_vector(0, n-1, err);			      forwardError(*err, __LINE__,);
   h = sm2_vector(0, n-1, err);			      forwardError(*err, __LINE__,);
   xi = sm2_vector(0, n-1, err);    		      forwardError(*err, __LINE__,);

   fp = post_cg(config, p, pb, quiet, err);           forwardError(*err, __LINE__,);

   dposterior_log_pdf(config, p, xi, err);            forwardError(*err, __LINE__,);

   for (j=0;j<n;j++) {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j];
   }
   for (its=1;its<=ITMAX;its++) {
      *iter=its;
      linmin(config,p,xi,fret,pb,quiet,err);                forwardError(*err, __LINE__,);

      if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+ZEPS)) {
	 goto end;
      }
      fp= *fret;
      dposterior_log_pdf(config,p,xi,err);	      forwardError(*err, __LINE__,);
      dgg=gg=0.0;
      for (j=0;j<n-1;j++) {
	 gg += g[j]*g[j];
	 dgg += (xi[j]+g[j])*xi[j];
      }
      if (gg == 0.0) {
	 goto end;
      }
      gam=dgg/gg;
      for (j=0;j<n-1;j++) {
	 g[j] = -xi[j];
	 xi[j]=h[j]=g[j]+gam*h[j];
      }
   }
   *err = addError(mcmc_tooManySteps, "too many steps", *err, __LINE__);

 end:
   sm2_free_vector(xi, 0, n-1);
   sm2_free_vector(h, 0, n-1);
   sm2_free_vector(g, 0, n-1);
   return;
}
#undef ITMAX
#undef ZEPS

/* ============================================================ *
 * Looks for the maximum-log-posterior point (stored in p) and  *
 * returns min(-logP).						*
 * ============================================================ */

double max_logP(config_base *config, double *p, const parabox *pb, double ftol,
		int quiet, error **err)
{
   double val;
   int iter;

   /* Minimization, cg-method (from NR) */
   frprmn(config, p, ftol, &iter, &val, pb, quiet, err);
   forwardError(*err, __LINE__, 0);
   return val;
}

void set_start_parameter(double *pstate, start_t start, const double *fid, config_base config, gsl_rng *rng,
			 FILE *FLOG, error **err)
{
   int i;

   if (FLOG) fprintf(FLOG, "p_ini =");

   for (i=0; i<config.npar; i++) {
      switch (start) {
	 case start_ran :
	    pstate[i] = gsl_ran_flat(rng, config.min[i], config.max[i]);
	    break;
	 case start_fid :
	    testErrorRet(fid==NULL, mcmc_undef, "Fiducial parameter ('fid' in config file) not set", *err, __LINE__,);
	    pstate[i] = fid[i];
	    break;
	 case start_min :
	    pstate[i] = config.min[i];
	    break;
	 case start_max :
	    pstate[i] = config.max[i];
	    break;
	 default :
	    *err = addErrorVA(mcmc_undef, "Unknown or invalid start point type %d(%s)",
			      *err, __LINE__, start, sstart_t(start));
	    return;
      }
      if (FLOG) fprintf(FLOG, " %g", pstate[i]);
   }
   if (FLOG) {
      fprintf(FLOG, "\n");
      fflush(FLOG);
   }
}

/* Returns 1 if maxP is maximum-posterior around parameter pstate */
int test_maximum(const double *pstate, const config_base config, double ffh, double maxP, error **err)
{
   int i, j, pm[2] = {+1, -1}, ok;
   double val, dx, *x;

   x = malloc_err(sizeof(double)*config.npar, err);   forwardError(*err, __LINE__, 1);
   for (i=0; i<config.npar; i++) {
      x[i] = pstate[i];
   }

   ok = 1;
   for (i=0; i<config.npar; i++) {
      for (j=0; j<2; j++) {
	 dx   = ffh*(config.max[i]-config.min[i]);
	 x[i] = pstate[i] + pm[j]*dx;

	 val = posterior_log_pdf_common_void((void*)(&config), x, err);
	 forwardError(*err, __LINE__, 1);
	 fprintf(stderr, "logPost for p%d: %s%s%g = %g", i, spar_t(config.par[i]), j==0?"+":"-", dx, val);
	 if (val>maxP) {
	    ok = 0;
	    fprintf(stderr, " (not ok)");
	 }
	 fprintf(stderr, "\n");
      }
      x[i] = pstate[i];
   }

   return ok;
}

/* ============================================================ *
 * Returns the log-likelihood values for the (starting) simplex *
 * vector of vectors p.						*
 * ============================================================ */
double *amoeba_starting_values(double **p, int npar, config_base *config, const parabox *pb,
			       int quiet, error **err)
{
   int i;
   double *y;

   y = malloc_err(sizeof(double)*(npar+1), err);  forwardError(*err, __LINE__, NULL);

   for (i=0; i<=npar; i++) {
      y[i] = post_cg(config, p[i], pb, quiet, err);
      forwardError(*err, __LINE__, NULL);
   }

   if (! quiet) fprintf(stderr, "done\n");
   return y;
}

/* ============================================================ *
 * The amoeba (simplex) minimisation method from Numerical      *
 * Recipes.							*
 * ============================================================ */
double amotry(double **p, double y[], double psum[], int ndim,
	      int ihi, double fac, config_base *config, const parabox *pb, int quiet, error **err)
{
   int j;
   double fac1,fac2,ytry,*ptry;

   ptry = malloc_err(sizeof(double)*ndim, err); forwardError(*err, __LINE__, 0.0);
   fac1 = (1.0-fac)/ndim;
   fac2 = fac1-fac;
   for (j=0;j<ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
   //ytry=(*funk)(ptry);
   ytry = post_cg(config, ptry, pb, quiet, err);
   forwardError(*err, __LINE__, 0.0);
   if (ytry < y[ihi]) {
      y[ihi]=ytry;
      for (j=0;j<ndim;j++) {
	 psum[j] += ptry[j]-p[ihi][j];
	 p[ihi][j]=ptry[j];
      }
   }
   free(ptry);
   return ytry;
}

#define TINY 1.0e-10
#define NMAX 5000
#define GET_PSUM \
  for (j=0;j<ndim;j++) {						\
     for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j];			\
     psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void amoeba(double **p, double y[], int ndim, double ftol, int *nfunk,
	    config_base *config, const parabox *pb, int quiet, error **err)
{
   int i,ihi,ilo,inhi,j,mpts=ndim+1;
   double rtol,sum,swap,ysave,ytry,*psum;

   psum = malloc_err(sizeof(double)*ndim, err);  forwardError(*err, __LINE__,);

   *nfunk=0;
   GET_PSUM
     for (;;) {
	ilo=0;
	ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
	for (i=0;i<mpts;i++) {
	   if (y[i] <= y[ilo]) ilo=i;
	   if (y[i] > y[ihi]) {
	      inhi=ihi;
	      ihi=i;
	   } else if (y[i] > y[inhi] && i != ihi) inhi=i;
	}
	rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
	if (rtol < ftol) {
	   SWAP(y[0],y[ilo])
	     for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i])
				    break;
	}
	testErrorRet(*nfunk >= NMAX, math_tooManySteps, "NMAX exceeded", *err, __LINE__,);
	*nfunk += 2;
	ytry = amotry(p,y,psum,ndim,ihi, -1.0, config, pb, quiet, err);
	forwardError(*err, __LINE__,);
	if (ytry <= y[ilo]) {
	   ytry = amotry(p,y,psum,ndim,ihi, 2.0, config, pb, quiet, err);
	   forwardError(*err, __LINE__,);
	} else if (ytry >= y[inhi]) {
	   ysave = y[ihi];
	   ytry  = amotry(p,y,psum,ndim,ihi, 0.5, config, pb, quiet, err);
	   forwardError(*err, __LINE__,);
	   if (ytry >= ysave) {
	      for (i=0;i<mpts;i++) {
		 if (i != ilo) {
		    for (j=0;j<ndim;j++)
		      p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
		    y[i] = post_cg(config, psum, pb, quiet, err);
		    forwardError(*err, __LINE__,);
		 }
	      }
	      *nfunk += ndim;
	      GET_PSUM
		}
	} else --(*nfunk);
     }

   free(psum);
}
#undef SWAP
#undef GET_PSUM
#undef NMAX
#undef TINY

/* ======================= -*- c -*- ========================== *
 * likeli.i                                                     *
 * MK 12/2006                                                   *
 * ============================================================ */


/* write, format="%s\n", "Usage:";
   write, format="   %s\n", "imin = read_chin(\"chi2\", n);";
   write, format="   %s\n", "marginn, chinall, k1, k2, which_chi2, which_prior;";
   write, format="   %s\n", "plot_chi2n, chi2, k1, k2, minpar=chinall(1:n-2,imin);";
   write, format="   %s\n", "// marginalize even if chi2 depends on only ";
   write, format="   %s\n", "// two parameters, to get the 2d-array chi2.";
*/

   
#include "stuff.i"
#include "spline.i"

ludo = 0;
jonben = 1;

func read_chin(name, n, which_chi2=, prior=, quiet=)
{
   /* DOCUMENT read_chin(name, n, which_chi2=, noprior=, quiet)
        Reads a chi^2 file with format
        p1 p2 p3 ... pn chi^2 chiz^2.
    */
   
   if (is_void(name)) error, "name of chi file?";
   if (is_void(n)) error, "number of parameters n?";
   if (is_void(which_chi2)) which_chi2 = 1;
   if (is_void(prior)) prior = 1;
   if (is_void(quiet)) quiet = 0;

   if (which_chi2!=1 && which_chi2!=2)
     error, "wrong which_chi2 (must be 1 or 2)";

   if (prior!=0 && prior!=1) error, "wrong prior (must be 0 or 1)\n";
   
   extern chinall, par, nel;

   N = numberoflines(name);
   if (N==-1) {
      mess = swrite(format="File %s not found", name);
      error, mess;
   }
   chinall = array(double, n+1+prior, N);
   nel = array(int, n);

   f = open(name);
   read, f, chinall;
   close, f;

   /* number of grid points for each dimension */
   for (i=n,p=1; i>=2; i--) {
      nel(i) = where(chinall(i,)==chinall(i,1))(1+p)/p;
      if (i==n) nel(i) -= 1;
      p *= nel(i);
   }
   nel(1) = N/p;

   par = array(double, n, max(nel));
   for (i=n,p=1; i>=1; i--) {
      par(i,1:nel(i)) = chinall(i,::p)(1:nel(i));
      p *= nel(i);
   }

   if (quiet==0) {
      write, format="%s", "nel = ";
      nel;
   }

   imin = where(min(chinall(n+which_chi2,))==chinall(n+which_chi2,));

   if (is_scalar(imin) && quiet==0) {
      write, format="minimum = %.2f for p=", chinall(n+which_chi2,imin);
      write, format="%.3f ", chinall(1:n,imin); writeNL;
   }
   return imin;
}

func get_chin(chinall, iel, which_chi2)
{
   /* DOCUMENT get_chin(chinall, iel, which_chi)
        Returns the chi^2-value corresponding to the index array iel.
        which_chi2 = 1,2.
    */

   if ((n=dimsof(iel)(2))!=dimsof(nel)(2))
     error, "wrong number of dimensions.";

   if (which_chi2!=1 && which_chi2!=2)
     error, "wrong which_chi2 (must be 1 or 2)";

   for (i=n, myi=0, p=1; i>=1; i--) {
      if (iel(i)>nel(i)) error, "index overreach";

      myi += (iel(i)-1)*p;
      p   *= nel(i);
   }

   myi++;

   return chinall(n+which_chi2, myi);
}

func next(&iel, k1, k2)
{
   /* DOCUMENT next(&iel, k1, k2)
        Increases the n-dim index array iel lexically by one. Returns
        the dimension furthest left that was changed, 0 if iel=nel.
        If k1, k2 are given those dimensions are unchanged.
    */

   if ((n=dimsof(iel)(2))!=dimsof(nel)(2))
     error, "wrong number of dimensions.";

   for (i=n,res=0; i>=1; i--) {
      if (!is_void(k1) && k1==i) continue;
      if (!is_void(k2) && k2==i) continue;
      
      if (iel(i)<nel(i)) {
         iel(i)++;
         res = i;
         break;
      } else {
         iel(i) = 1;
         continue;
      }
   }
   return res;
}

func fix_1d(chinall, k, v)
{
   /* DOCUMENT fix_1d(chinall, k, v)
        Returns chi^2 on the subspace where the k-th parameter equals v.
    */

   if (k<1 || k>dimsof(nel)(2)) error, "wrong column number";

   Npar = dimsof(nel)(2);

   whereisv = where2(chinall(k,)==v);
   if (!is_array(whereisv)) {
      //error, "parameter value not taken";
      write, format="%s\n", "interpolate to parameter value";
      i0 = where(chinall(k,)<v)(0);
      v0 = chinall(k,i0);
      v1 = chinall(k,i0+1);
      whereisv0 = where2(chinall(k,)==v0);
      whereisv1 = where2(chinall(k,)==v1);
      r  = (v-v0)/(v1-v0);
      r;
      chifix = (1-r)*chinall(where(indgen(Npar+2)!=k), whereisv0(1,)) + r*chinall(where(indgen(Npar+2)!=k), whereisv1(1,));
   } else {
      chifix   = chinall(where(indgen(Npar+2)!=k), whereisv(1,));
   }

   par = par(where(indgen(Npar)!=k),);
   nel = nel(where(indgen(Npar)!=k));
   
   return chifix;
}

func marginn(chinall, k1, k2, which_chi2, which_prior)
{
   /* DOCUMENT marginn(chinall, k1, k2, which_chi2, which_prior)
        Marginalises the n-dimensional chinall, with keeping the
        dimensions k1 and k2.
        which_chi2 = 1,2.
        which_prior = 0,1,2 (0: no prior).
    */

   big    = 1.0e10;
   thresh = big;

   if (which_chi2!=1 && which_chi2!=2) error, "wrong which_chi2 (must be 1 or 2)";
   if (which_prior<0 || which_prior>2) error, "wrong which_prior (must be 0, 1 or 2)";
   if (which_chi2==which_prior) error, "which_chi2=which_prior!";

   n = dimsof(nel)(2);
   extern chi2;
   chi2 = array(double, nel(k1), nel(k2));
   delta = par(,dif)(,1);
   
   for (i=1; i<=nel(k1); i++) {
      for (j=1; j<=nel(k2); j++) {

         iel = array(1, n);
         iel(k1) = i;
         iel(k2) = j;
         L   = 0.0;
         res = 1;
         
         do {

            x = -0.5*get_chin(chinall, iel, which_chi2);
            if (x>-thresh) expx = exp(x);
            else expx = 0.0;

            if (which_prior!=0) {
               prior = -0.5*get_chin(chinall, iel, which_prior);
               if (prior>-thresh) expx *= exp(prior);
            }

            L += expx;

            res = next(iel, k1, k2);

         } while (res>0);

         chi2(i,j) = L;
         for (k=1; k<=n; k++) {
	    /* Bug fixed (j2->k3) (21/09/09) */
            if (k==k1 || k==k2) continue;
            chi2(i,j) *= delta(k);
         }
         if (chi2(i,j)>0) chi2(i,j) = -2.0*log(chi2(i,j));
         else chi2(i,j) = big;
         
      }
   }
}

func marginn_1d(chinall, k, which_chi2, which_prior)
{
   /* DOCUMENT marginn_1d(chinall, k, which_chi2, which_prior)
        SEE ALSO marginn.
    */

   big    = 1.0e10;
   thresh = big;

   if (which_chi2!=1 && which_chi2!=2) error, "wrong which_chi2 (must be 1 or 2)";
   if (which_prior<0 || which_prior>2) error, "wrong which_prior (must be 0, 1 or 2)";
   if (which_chi2==which_prior) error, "which_chi2=which_prior!";

   n = dimsof(nel)(2);
   extern chi2;
   chi2 = array(double, nel(k));
   delta = par(,dif)(,1);

   for (i=1; i<=nel(k); i++) {
      iel = array(1, n);
      iel(k) = i;
      L      = 0.0;
      res    = 1;

      do {

         if (which_prior!=0) {
            prior = -0.5*get_chin(chinall, iel, which_prior);
            if (prior>-thresh) expp = exp(prior);
            else expp = 0.0;
         } else expp = 1.0;

         x = -0.5*get_chin(chinall, iel, which_chi2);
         if (x>-thresh) expx = exp(x);
         else expx = 0.0;

         L += expx*expp;

      afterL:

         res = next(iel, k);

      } while (res>0);

      chi2(i) = L;
      for (j=1; j<=n; j++) {
         if (j==k) continue;
         chi2(i) *= delta(j);
      }
      if (chi2(i)>0) chi2(i) = -2.0*log(chi2(i));
      else chi2(i) = big;

   }
      
}

func interpol2d(a, x, y, newx, newy)
{
   /* DOCUMENT interpol2d(a, x, y, newx, newy)
        Interpolates a given on (x,y) onto the new 2d array
        (newx,newy).
    */

   dx = x(2)-x(1);
   dy = y(2)-x(1);

   Nx    = dimsof(x)(2);
   Ny    = dimsof(y)(2);
   Nnewx = dimsof(newx)(2);
   Nnewy = dimsof(newy)(2);
   newa  = array(double, Nnewx, Nnewy);

   for (i=1; i<=Nnewx; i++) {
      tx = (newx(i)-x(1))/dx;
      ix = int(floor(tx));
      dix = tx-ix;
      ix++;
      if (ix<1 || ix>Nx-1) continue; //error, "overflow x";
      for (j=1; j<=Nnewy; j++) {
         ty = (newy(j)-y(1))/dy;
         jy = int(floor(ty));
         djy = ty-jy;
         jy++;
         if (jy<1 || jy>Ny-1) continue; //error, "overflow y";

         //[dix, diy];
         newa(i,j) = (1.0-dix)*(1.0-djy)*a(ix,jy) \
           + (1.0-dix)*djy*a(ix,jy+1) \
           + dix*(1.0-djy)*a(ix+1,jy) \
           + dix*djy*a(ix+1,jy+1);

      }
   }
   return newa;
}

func add_key(key, key_str, Col, width, type, height)
{
   if (!is_void(key) && key!=0) {
      l  = limits();
      dx  = l(2)-l(1);
      dy  = l(4)-l(3);

      /* Lines */
      fxl = 0.75*dx;
      dxl = 0.15*dx;
      fyl = 0.9*dy;
      dyl = 0.075*dy;

      xl = l(1) + fxl;
      x  = [xl, xl+dxl];
      yl = l(3) + fyl;
      y  = [yl, yl]-(key-1)*dyl;
      plg, y, x, marks=0, color=Col, width=width, type=type;

      /* Label */
      dxk = 0.05*dx;
      if (is_void(height)) { height = 24; }
      plt, key_str, xl-dxk, yl-(key-1)*dyl, justify="RH", tosys=1, color=Col, height=height / 1.2;
   }
}

func plot_chi2n(chi2, k1, k2, alllevels=, minpar=, nice=, do_fma=,
                type=, chi2thres=, Omegam0=, fitpar=, Col=, shade=, width=, quiet=, levels1=, key=, key_str=, height=)
{
/* DOCUMENT plot_chi2n(chi2, k1, k2, alllevels=, minpar=, nice=, do_fma=, type=,
      chi2thres, Omegam0=, fitpar=, Col=, shade=, width=, quiet=, levels1=, key=)
      Contour plot of chi2. alllevels=0: 1,2,3sigma, 1: all kinds of levels,
      2: 1sigma using color Col(,1), 3: 1,2sigma using colors Col(,1), Col(,2).
      Own levels (levels1) can be specified (alllevels should then be 0).
 */

   if (is_void(alllevels)) alllevels = 0;
   if (is_void(nice)) nice = 1;
   if (is_void(do_fma)) do_fma = 1;
   if (is_void(type)) type = "solid";
   if (is_void(chi2thres)) chi2thres = 0;
   if (is_void(Omegam0)) Omegam0 = 0.3;
   if (is_void(width)) width = 4;
   if (is_void(quiet)) quiet = 0;

   
   extern gx, gy, m, par1, par2;
   col = [Red, Darkgreen, Blue, Orange, Violett, Cyan];

   if (is_void(Col)) Col = col;

   if (is_void(k1)) k1 = 1;
   if (is_void(k2)) k2 = 2;

   n1 = nel(k1);
   n2 = nel(k2);

   par1 = par(k1, 1:n1);
   par2 = par(k2, 1:n2);

   gx = par1(,-:1:n2);
   gy = par2(-:1:n1,);
     
   /* Frequentist chi^2 levels for 2 dof (see NR p.697) */
   if (is_void(levels1)) {
      levels1 = [2.3, 6.17, 11.8];            /* 1, 2, 3 sigma       */
   }
   levels2 = [4.61, 9.21, 18.4];              /* 1.5, 2.5, 3.5 sigma */
   levels3 = [20, 50, 100];
   levels4 = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];

   logxy, 0, 0;
   if (do_fma) fma;

   m = min(chi2);
   w = where2(chi2==m);
   if (quiet==0) {
      write, format="minimum of marginalized chi^2: %.3f for (%.3f, %.3f)\n", m, par(k1,w(1)), par(k2,w(2));
   }

   smooth = 1;
   if (nice==1) {
      if (!is_void(shade)) {
         set_palette_shade, shade, Npal=200;
	 /* TODO: Factor in exponential depending on histogram size or smoothing factor */
         e = exp(-chi2*3);
         f = 1.0001; // f = 1.001;
         plfc, e, gy, gx, levs=span(min(e), max(e)/f, 100);
      }
      if (alllevels<2) {
         for (s=1; s<=dimsof(levels1)(2); s++) {
            plc, chi2-m, gy, gx, levs=levels1(s), smooth=smooth, marks=0, width=width, color=Col(,s), type=type;
         }
      }
      if (alllevels>=2) {
         plc, chi2-m, gy, gx, levs=levels1(1), smooth=smooth, marks=0, width=width, color=Col(,1), type=type;
      }
      if (alllevels>=3) {
         plc, chi2-m, gy, gx, levs=levels1(2), smooth=smooth, marks=0, width=width, color=Col(,2), type=type;
      }
      if (alllevels>=4) {
         plc, chi2-m, gy, gx, levs=levels1(3), smooth=smooth, marks=0, width=width, color=Col(,3), type=type;
      }
   } else {
      plc, chi2-m, gy, gx, levs=levels1, smooth=smooth, marks=0, width=width, type=type;
   }

   if (chi2thres>0) {
      w = where2(chi2-m<chi2thres);
      //plmk, par2(w(2,)), par1(w(1,)), marker=4, width=10, msize=0.25;
      l = limits();
      x = span(l(1), l(2), 50);
      plg, fitpar(3)*(x/Omegam0)^-fitpar(1), x, marks=0, width=width;
   }
   
   /* palette, "stern.gp"; */
   /* plfc, chi2-m, gy, gx, levs=levels1; */
   if (alllevels==1) {
      plc, chi2-m, gy, gx, levs=levels2, smooth=smooth, marks=0, type="dot";
      plc, chi2-m, gy, gx, levs=levels3, smooth=smooth, marks=0, type="dash";
      plc, chi2-m, gy, gx, levs=levels4, smooth=smooth, marks=0, type="dashdot";
   }

   /* plot grid */
   /* plm, gy, gx, type="dot"; */

   if (!is_void(minpar) && nice==0) {
      /* plmk, minpar(k2), minpar(k1), marker=4, color=Red, width=10, msize=0.5; */
      plmk, par(k2,w(2)), par(k1,w(1)), marker=4, color=Red, width=10, msize=0.5;
   }

   add_key, key, key_str, Col(,1), width, type, height;
}

func plot_chi2_1d(chi2_1d, k, &meansig, &pmax, outputfile=, nice=, do_fma=, do_mean=, do_sigma=,
                  color=, Type=, mcmc=, text=, width=, norm=, key=, key_str=, height=)
{
   /* DOCUMENT plot_chi2_1d(chi2_1d, k, &meansig, &pmax, outputfile=, nice=, do_fma=,
                            do_mean=, do_sigma=, type=, mcmc=, text=, width=, norm=,
			    key=, key_str=)
        
    */

   if (is_void(k)) k = 1;
   if (is_void(nice)) nice = 0;
   if (is_void(do_fma)) do_fma = 0;
   if (is_void(do_mean)) do_mean = 1;
   if (is_void(do_sigma)) do_sigma = [1];
   if (is_scalar(do_sigma)) do_sigma = [do_sigma];
   if (is_void(color)) color = Black;
   if (is_void(mcmc)) mcmc = 0;
   if (is_void(text)) text = "";
   if (is_void(norm)) norm = 0; /* Maximum normalisation */


   logxy, 0, 0;
   if (do_fma) fma;

   N = dimsof(chi2_1d)(2);

   if (mcmc==0) L = exp(-0.5*(chi2_1d-min(chi2_1d)));
   else if (mcmc==1) {
      if (norm==0) { L = chi2_1d/max(chi2_1d); }
      else if (norm==1) {
	 delta = par(k,2)-par(k,1);
	 L     = chi2_1d/(sum(chi2_1d)*delta);
      }
      else {
	 error, "Undefined normalisation mode (must be 0 for maximum or 1 for integral norm)";
      }
   }
   
   if (nice==0) {
      /* Data points */
      plmk, L, par(k, 1:nel(k)), marker=4, msize=0.3, width=10;
   }

   /* Spline */
   Nsp    = 5;
   extern parsp, Lsp;
   parsp  = span(par(k,1), par(k,N), nel(k)*Nsp);
   if (mcmc==0) {
      /* TODO: different max treatment, as for mcmc */
      chi2sp = spline(chi2_1d, par(k, 1:nel(k)), parsp);
      Lsp = exp(-0.5*(chi2sp-min(chi2sp)));
   } else if (mcmc==1) {
      /* Tension goes between 0 (= cubic spline) to infinity (linear interp) */
      //Lsp = tspline(9, chi2_1d/max(chi2_1d), par(k, 1:nel(k)), parsp);
      Lsp = tspline(9, L, par(k, 1:nel(k)), parsp);
   }

   /* Maximum parameter */
   wmax = where(Lsp==max(Lsp))(1);
   pmax  = parsp(wmax);

   /* Plot posterior/likelihood */
   plg, Lsp, parsp, marks=0, color=color, type=Type, width=width;

   if (do_mean) {
      /* mean */
      mean     = mean_1d(L, k, Nspline=0);
      mean_sp  = mean_1d(Lsp, k, Nspline=Nsp);
      //plg, [0, 1], [mean, mean], marks=0, color=color;
      plg, [0, interp(Lsp, parsp, mean_sp)], [mean_sp, mean_sp], marks=0, type=2,
            color=color;
      if (!is_void(meansig)) meansig(1) = mean_sp;
   }

   if (!is_void(do_sigma)) {
      /* n-sigma region */
      for (s=1; s<dimsof(do_sigma)(2); s++) {
         wsigma = do_sigma(s);
         if (wsigma==0) { continue; }
         sigma_pm = sigma_1d(L, k, Nspline=0, wsigma=wsigma);
         sigma_pm_sp = sigma_1d(Lsp, k, Nspline=Nsp, wsigma=wsigma);
         plg, [0, interp(Lsp, parsp, sigma_pm_sp(1))], [sigma_pm_sp(1,), sigma_pm_sp(1,)],
           marks=0, type=3, color=color;
         plg, [0, interp(Lsp, parsp, sigma_pm_sp(2))], [sigma_pm_sp(2,), sigma_pm_sp(2,)],
           marks=0, type=3, color=color;
         if (!is_void(meansig)) meansig(2:3) = sigma_pm_sp;
         if (wsigma==1) { one_sigma_pm_sp = sigma_pm_sp; }
      }
   }

   //limits;

   if (!is_void(outputfile) && do_mean==1 && !is_void(one_sigma_pm_sp)) {
      f = open(outputfile, "w");
      write, f, format="# mean pm %dsigma (Nspline=%d) # %s\n", wsigma, Nsp, text;
      write, f, format="%.4f + %.4f - %.4f\n", mean_sp,
        (mean_sp-one_sigma_pm_sp(1)), (one_sigma_pm_sp(2)-mean_sp);
   }

   add_key, key, key_str, color, width, height;
}

func mean_1d(L, k, Nspline=)
{
   /* DOCUMENT mean_1d(L, k, Nspline)
        Returns the mean = int(p*L)/int(L)).
        New call: L instead of chi2!
    */

   if (is_void(Nspline)) Nspline = 0;

   if (Nspline==0) {
      p = par(k,1:dimsof(chi2_1d)(2));
   } else {
      p  = span(par(k,1), par(k,0), nel(k)*Nspline);
   }
   
   /* Integrate over p*L from first to last data point */
   return integ(p*L, p, p(0))/integ(L, p, p(0));
}

func sigma_1d(L, k, Nspline=, wsigma=)
{
   /* DOCUMENT sigma_1d(chi2_1d, k, Nspline=0)
        Returns [+1sigma, -1sigma].
        New call: L instead of chi2!
    */

   if (is_void(Nspline)) Nspline = 0;
   if (is_void(wsigma)) wsigma = 1;
   if (wsigma!=1 && wsigma!=2 && wsigma!=3) error, "which sigma?";

   sig = [0.683, 0.954, 0.9973];
   /* 68.3%, ... */
   
   if (Nspline==0) {
      p = par(k,1:dimsof(L)(2));
   } else {
      p  = span(par(k,1), par(k,0), nel(k)*Nspline);
   }

   sigma_pm = array(0.0, 2);
   N = dimsof(L)(2);
   for (i=1; i<=N; i++) {
      x_m = integ(L, p, p(i))/integ(L, p, p(0));
      if (x_m>=(1.0-sig(wsigma))/2.0) {
         sigma_pm(1) = p(i);
         break;
      }
   }
   LL = L(::-1);
   pp = p(::-1);
   for (i=1; i<=N; i++) {
      x_p = integ(LL, pp, pp(i))/integ(LL, pp, pp(0));
      if (x_p>=(1.0-sig(wsigma))/2.0) {
         sigma_pm(2) = pp(i);
         break;
      }
   }

   return sigma_pm;
}

func log_power_law(logx, p)
{
   /* DOCUMENT log_power_law(logx, p)
        y = p(1)*x^p(2). Returns logy.
    */
   return p(1) + p(2)*logx;
}

func fit_power_law(chi2, k1, k2, chi2thres, Omegam0=, outname=, mcmc=)
{
   if (is_void(Omegam0)) Omegam0 = 0.3;
   if (!is_void(outname) && outname!="") {
      outfile = open(outname, "w");
   }
   if (is_void(mcmc)) mcmc = 0;

   if (is_void(k1)) k1 = 1;
   if (is_void(k2)) k2 = 2;

   par1 = par(k1, 1:nel(k1));
   par2 = par(k2, 1:nel(k2));

   if (mcmc==1) {
      /* Density -> chi^2 */
      extern ctmp;
      ctmp = chi2;
      w = where(chi2);
      ctmp(w) = -2.0*log(chi2(w));
      w = where(chi2==0);
      ctmp(w) = 1.0e10;
      chi2 = ctmp;
   }

   m = min(chi2);
   w = where2(chi2-m<chi2thres);
   write, format="Fitting power law from %d data points (chi2thres=%.2f)\n", dimsof(w)(3), chi2thres;
   if (!is_void(outfile)) write, outfile, format="Fitting power law from %d data points\n", dimsof(w)(3);
   
   extern result, l1, l2, p, alpha, deltaalpha, s8, deltas8;
   l1 = log(par1(w(1,)));
   l2 = log(par2(w(2,)));
   /* p(1) = s8, p(2) = log(index) */
   p = [1.0, 1.0];
   weight = array(1.0, dimsof(l2)(2));
   for (i=1; i<=dimsof(l2)(2); i++) {
      c2 = chi2(w(1,i), w(2,i))-m;
      weight(i) = exp(-0.5*c2);
   }

   /* yorick:yutils have to be installed */
   include, "yutils/lmfit.i";
   result = lmfit(log_power_law, l1, p, l2, weight, correl=1, stdev=1);

   alpha      = -p(2);
   deltaalpha = (*result.stdev)(2);
   s8         = exp(p(1))*Omegam0^p(2);
   deltas8    = exp(p(1))*(*result.stdev)(1)*Omegam0^p(2);

   write, format="Best fit: sigma_8 ( Omegam/%.2f )^( %.3f +- %.3f ) = %.4f +- %.4f\n",
     Omegam0, alpha, deltaalpha, s8, deltas8;
   if (!is_void(outfile)) write, outfile, format="Best fit: sigma_8 ( Omegam/%.2f )^( %.3f +- %.3f ) = %.4f +- %.4f\n",
     Omegam0, alpha, deltaalpha, s8, deltas8;
   return [alpha, deltaalpha, s8, deltas8];
}
   
func labels(t, x, y)
{
   pltitle, t;
   xytitles, x, y;
}

func plot_chi2z(nofz)
{
   /* DOCUMENT plot_chi2z(void)
        Plots chi^2 of the redshift parameters, for chi^2 file with
        alpha_p, beta_p, z0.
    */


   if (nofz==ludo) labels = ["!a_p", "!b_p", "z_0"];
   else if (nofz==jonben || nofz==ymmk) labels = ["a", "b", "c"];
   else error, "nofz=? [ludo, jonben, ymmk]";

   names = ["ab", "ac", "bc"];

   imin = read_chin("chi2", 3, which_chi2=2);

   alllevels = 0;
   Col = [Red, Blue, Orange];

   /* ab */
   marginn, chinall, 1, 2, 2, 0;
   plot_chi2n, chi2, 1, 2, minpar=chinall(1:3,imin), nice=1, alllevels=alllevels, Col=Col;
   xytitles, labels(1), labels(2);
   names(1); limits; limits();
   eps, names(1);

   /* ac */
   marginn, chinall, 1, 3, 2, 0;
   plot_chi2n, chi2, 1, 3, minpar=chinall(1:3,imin), nice=1, alllevels=alllevels, Col=Col;
   xytitles, labels(1), labels(3);
   names(2); limits; limits();
   eps, names(2);

   /* bc */
   marginn, chinall, 2, 3, 2, 0;
   plot_chi2n, chi2, 2, 3, minpar=chinall(1:3,imin), nice=1, alllevels=alllevels, Col=Col;
   xytitles, labels(2), labels(3);
   names(3); limits; limits();
   eps, names(3);
}

func out_chi2(chi2, chi2name, k1, k2)
{
   /* DOCUMENT out_chi2(chi2, chi2name, k1, k2)
        Writes the 2d-array chi2 to the file "chi2name".
    */

   par1 = par(k1, 1:nel(k1));
   par2 = par(k2, 1:nel(k2));
   gx = par1(,-:1:nel(k2));
   gy = par2(-:1:nel(k1),);

   f = open(chi2name, "w");
   write, f, gx, gy, chi2-min(chi2);
   close, f;
}

func read_all_proposals(niter, dir=, name=)
{
   /* DOCUMENT read_all_proposals(niter, dir=, name=)
        Reads mix_mvdens files name_i for i=0..niter and returns array of struct mix_mvdens.
    */

   if (is_void(niter)) error, "niter=?";
   
   if (is_void(dir)) dir = ".";
   if (is_void(name)) name = "proposal";

   mix_mvd = array(mix_mvdens, niter+1);    /* 1 .. niter+1 */
   for (iter=1; iter<=niter+1; iter++) {
      pname = swrite(format="%s/%s_%d", dir, name, iter-1);
      //write, format="reading %s\n", pname;

      mix_mvd(iter) = mix_mvdens_read(pname);
   }

   return mix_mvd;
}

func plot_mix_mvdens_mean(mix_mvd, xy, msize=, color=, connect=)
{
   /* DOCUMENT plot_mix_mvdens(mix_mvd, msize=, color=, connect=)
        Plots array of mix_mvdens as (connected) circles for their mean.
      SEE ALSO from read_all_proposals
    */

   min_weight = 1.0e-5;
  
   nmix  = dimsof(mix_mvd)(2);
   ncomp = mix_mvd.ncomp(1);

   if (is_void(msize)) msize = 1;

   c_frac = limits()/5;
   msize  = c_frac(2)-c_frac(1);
   msizey = c_frac(4)-c_frac(3);

   if (is_void(color)) {
      color = color_gnu;
   }
   if (is_void(connect)) connect = 0;
   else connect = 1;
   
   means = array(double, nmix, 2);
   flag  = array(int, ncomp, nmix);

   for (c=1; c<=ncomp; c++) {

      for (m=1; m<=nmix; m++) {

         mvd       = *(mix_mvd(m)).mvd;
         means(m,) = (*(mvd(c).mean))(xy);
         weight    = (*(mix_mvd(m).weight))(c);
         if (weight<min_weight) {
            flag(c,m) = 0;
            continue;
         }
         flag(c,m) = 1;

	 //plmk, means(m,2), means(m,1), marker=4, msize=weight, width=11, color=color(c);
	 plot_circle, means(m,[1,2]), size=msize*weight, sizey=msizey*weight, color=color(c%Ncolor_gnu),
	   width=(m==1?4:1);
      }

      wh = where(flag(c,));
      if (is_array(wh) && connect==1) plg, means(wh,2), means(wh,1), marks=0, color=color(c);
   }
}

func plot_mix_mvdens_covar(mix_mvd, xy, color=, plev=)
{
    /* DOCUMENT plot_mix_mvdens_covar(mix_mvd, xy, plev=)
          Plots 1,2,3sigma-contours of a mix_mvdens (each component)
     */

    ncomp = mix_mvd.ncomp;

    if (is_void(color)) {
      color = -indgen(4:4+ncomp);
      color(1) = -3;  /* black */
   }

   if (is_void(plev)) plev = 0;

   extern Z;
   for (c=1; c<=ncomp; c++) {
      weight = (*(mix_mvd.weight))(c);
      if (weight==0) continue;
      mvd = (*mix_mvd.mvd)(c);

      /* TODO: Is this really the right thing to do?? */
      /* Why are the eigenvalues negative sometimes? */
      mean = (*mvd.mean)(xy);
      var  = (*mvd.var)(xy, xy);
      for (s=1; s>=1; s--) {
	 sigma = sqrt(-2*log(1.0-erf(s/sqrt(2))));
	 Z = ellipse(mean, var, sigma=sigma);
	 if (!is_void(Z)) draw_ellipse, Z, type=2, color=color(c);
      }
      if (plev==1) {
	 str = swrite(format="(%.4f,%.4f)", eigenvector(1,1), eigenvector(2,1));
	 l = limits();
	 x = l(2)-0.05*(l(2)-l(1));
	 y = l(4)-0.05*(l(4)-l(3));
	 plt, str, x, y, justify="RT", tosys=1;
      }
   }

}

func plot_prop_weights_iter(mix_mvd, color=)
{
   /* DOCUMENT plot_prop_weight_iter(mix_mvd)
        Plots the weights of the proposal components as function of the iteration.
    */

   nmix  = dimsof(mix_mvd)(2);
   ncomp = mix_mvd.ncomp(1);

   if (is_void(color)) {
      color = -indgen(4:4+ncomp);
      color(1) = -3;  /* black */
   }

   weights = array(double, nmix, ncomp);
   for (m=1; m<=nmix; m++) {
      weights(m,) = *(mix_mvd.weight(m));
   }

   iter = indgen(nmix)-1;
   for (c=1; c<=ncomp; c++) {
      plg, weights(,c), iter, marks=0, color=color(c);
   }
   
   range, -0.05, 1.05;
   xytitles, "iteration", "component weight";
}

func max_hist(file)
{
   /* DOCUMENT max_hist(file)
        Calculates the maximum of a histogram using
	cubic spline interpolation.
    */

   Nsp = 10;

   read_table, file, hist;
   N = dimsof(hist)(3);
   xp = span(hist(1,1), hist(1,N), N*Nsp);
   yp = spline(hist(2,), hist(1,), xp);
   w = where(yp==max(yp));
   return xp(w);
}

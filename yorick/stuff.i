/* ======================= -*- c -*- ========================== *
 * stuff.i                                                      *
 * MK 2004-2007                                                 *
 * Can be included in other yorick files.                       *
 * ============================================================ */


/* ============================================================ *
 * About vectors and matrices:                                  *
 *   - vectors are column vectors                               *
 *   - the inner most brackets of a matrix denotes its          *
 *     columns                                                  *
 *                                        | 1  2 |              *
 *      e.g. [[1,3],[2,4]] is the matrix  |      |              *
 *                                        | 3  4 |              *
 *                                                              *
 *   - indices of matrices are (Zeile, Spalte)                  *
 *   - matrix multiplication is a(,+)*b(+,)                     *
 * ============================================================ */


/* ============================================ *
 * Some constants.                              *
 * ============================================ */

extern twopi, arcmin;

twopi     = 2.*pi;
arcmin    = 2.90888e-4;


/* ============================================ *
 * Useful functions.                            *
 * ============================================ */

func numberoflines(name, &ncomment)
{
   /* DOCUMENT numberoflines(name, &ncomment)
         Returns the number of lines of file name.
    */

   f = open(name, "r", 1);
   if (!f) return -1;

   i = ncomment = array(long);
   while (s = rdline(f)) {
      if (strpart(s, 1:1)=="#") {
         ncomment++;
      }
      i++;
   }
   close, f;
   return i;
}

func numberofcolumns(name, all=)
{
   if (is_void(all)) all = 1;
   
   NCOLMAX = 1024;
   dummy = array(string, NCOLMAX);
   f = open(name, "r", 1);
   if (!f) return -1;

   c = -1;
   while (s = rdline(f)) {
     if (strpart(s, 1:1)=="#") continue;
     cc = sread(s, format="%s", dummy);
     if (c==-1) c = cc;
     if (c!=cc) error, "file has not constant number of columns";
     if (all==0) break;
   }
   if (c==NCOLMAX) error, "overflow";
   return c;
}

func pwd(void)
{
   /* DOCUMENT pwd(void)
        Returns the present working directory.
    */

   return get_cwd();
}

func win(number, style=)
{
   /* DOCUMENT win(number, style=)
         Pops up a boxed window or switches to window # number if it already exists.
    */

   if (is_void(style)) style = "boxed-bigaxes.gs";

   if (current_window()!=number)
     window, number, style=style;
   else window, number;

   logxy, 1, 1;
}

func copylimits(win)
{
   /* DOCUMENT copylimits(win)
        Copies the limits of window win to the present window.
    */

   cw = current_window();
   window, win;
   lim = limits();
   window, cw;
   limits, lim(1), lim(2), lim(3), lim(4);
}

func abssqr(v)
{
   for (s=0,i=1; i<=dimsof(v)(2); i++) {
      s+=v(i)^2;
   }
   return s;
}

/* ============================================================ *
 * Astronomical calculations.                                   *
 * ============================================================ */

func deg2rad(h, m, s, ra=)
{
   while (is_void(ra) && (ra!=0 && ra!=1)) {
      write, format="%s\n", "ra/dec? [1/0]";
      ra = array(int);
      read, ra;
   }

   rad = h + m/60.0 + s/3600.0;
   if (ra==1) rad *= 15.0;

   return rad;
}
  


/* ============================================================ *
 * Colors.							*
 * ============================================================ */

/* === Some self-defined color with speaking names */
/* r g b */
White   = [255,255,255];
Black   = [0,0,0];
Blue    = [0,0,255];
Red     = [255,0,0];
Green   = [0,255,0];
Cyan    = [0,255,255];
Yellow  = [255,255,0];
Magenta = [255,0,255];
Lightblue = [150,200,255];
Darkgreen = [50,170,50];
Darkgreen2 = [50,140,50];
Violett = [250,50,250];
Orange  = [255,150,0];
Darkred = [170,30,00];

Color   = [Black, Blue, Red, Orange, Lightblue, Darkgreen, Violett];
/* to be used as Color(,i) ! */

GreyBlue = [200,200,255];
GreyRed  = [255,150,150];
GreyGreen = [200,255,200];
GreyCyan = [200,255,255];
GreyOrange = [255,220,120];

/* === color_gnu are the gnuplot colors === */
/* black, red, green, blue, magenta, cyan, yellow */
color_gnu = [-3, -5, -6, -7, -9, -8, -10];
Ncolor_gnu = dimsof(color_gnu)(2);

func set_palette_shade(color, Npal=)
{
   /* DOCUMENT palette_shade(color, Npal=)
        Sets the palette to contain shades from color "color" to white.
    */

   if (is_void(Npal)) Npal = 100;

   Min = 50;

   if (color=="Blue") {
      red   = span(Min, 255, Npal);
      green = span(Min, 255, Npal);
      blue  = span(255, 255, Npal);
   } else if (color=="Red") {
      red   = span(255, 255, Npal);
      green = span(Min, 255, Npal);
      blue  = span(Min, 255, Npal);
   } else if (color=="Green") {
      red   = span(Min, 255, Npal);
      green = span(255, 255, Npal);
      blue  = span(Min, 255, Npal);
   } else if (color=="Darkgreen") {
      red   = span(Min, 255, Npal);
      green = span(170, 255, Npal);
      blue  = span(Min, 255, Npal);
   } else if (color=="Orange") {
      red   = span(255, 255, Npal);
      green = span(150, 255, Npal);
      blue  = span(0, 255, Npal);
   } else error("color (not yet) supported");

   palette, red, green, blue;
}

/* ============================================================ *
 * Matrix stuff.						*
 * ============================================================ */

func SVana(C, Nsingular)
{
   extern sigma, U, Vt;

   sigma = SVdec(C, U, Vt);

   if (Nsingular>dimsof(C)(2)) error, "Nsingular too large";

   U  = U(,1:Nsingular);
   Vt = Vt(1:Nsingular,);
   sigma = sigma(1:Nsingular);
   Sigma = make_diag(sigma);

   Cp = (U(,+)*Sigma(+,))(,+)*Vt(+,);
   Cp = Cp(1:Nsingular, 1:Nsingular);

   write, format="Nsingular = %d: ", Nsingular;
   if (cholesky(Cp)==0)  write, format="%s\n", "positive definite";
   else write, format="%s\n", "not positive definite";
   
   return Cp;
}

func indgen2d(n, &x1, &x2)
{
   /* DOCUMENT indgen2d(n)
        Generates and returns a 2xn matrix a which contains indices for a nxn matrix.
	x1 and x2 are optional.
    */

   a = array(int, 2, n^2);
   x = indgen(n);

   X1 = x()(,-:1:n);
   X2 = x()(-:1:n,);
   a(1,) = X1(where(X1));
   a(2,) = X2(where(X2));

   x1 = X1;
   x2 = X2;

   return a;
}

func matrix(a, n=)
{
   /* DOCUMENT matrix(a, n=)
      Creates a matrix from the n(1)*n(2)-dim vector a. Opposite of x(where(x))
      (for x!=0). Indices right??
    */

   if (is_void(n)) {
      n = array(int, 2);
      n(1) = n(2) = int(sqrt(dimsof(a)(2)));
   }
   
   A = array(double, n(1), n(2));

   for (k=1; k<=n(1); k++) {
      for (l=1; l<=n(2); l++) {
	 //j = k + (l-1)*n(1);
         j = (k-1)*n(2) + l;
	 A(k,l) = a(j);
      }
   }

   return A;
}

func diag_element(u, i)
{
   /* DOCUMENT diag_element(i, n)
        Returns the diagonal element (i,i,...,i) for the n-dimensional
        array u. Used by diag.
    */

   n = dimsof(u)(1);

   if (n==2) return u(i,i);
   if (n==3) return u(i,i,i);
   if (n==4) return u(i,i,i,i);
   if (n==5) return u(i,i,i,i,i);
   if (n==6) return u(i,i,i,i,i,i);
   error("wrong n in diag_element");
}

func diag(u)
{
   /* DOCUMENT diag(u)
         Returns the diagonal of u. u may have (in principle)
         arbitrary type and number of dimensions.
      SEE ALSO: make_diag, diag_matrix.
    */

   n = dimsof(u)(2);

   if (sizeof(diag_element(u,1))==4) d = array(float, n);
   else if (sizeof(diag_element(u,1))==8) d = array(double, n);
   else d = array(complex, n);

   for (i=1; i<=n; i++) d(i) = diag_element(u, i);

   return d;
}

func make_diag(x)
{
   /* DOCUMENT make_diag(x)
        Returns a diagonal matrix with the vector x on the diagonal.
        Thus, make_diag(diag(M)) make a diagonal matrix from M.
      SEE ALSO: diag, diag_matrix.
    */

   n = dimsof(x)(2);

   A = array(double, n, n);
   for (i=1; i<=n; i++) A(i,i) = x(i);
   return A;
}

func diag_matrix(A)
{
   /* DOCUMENT diag_matrix(A)
      Returns matrix A with off-diagonal set to zero.
      SEE ALSO: make_diag, diag
    */

   n = dimsof(A)(2);

   return unit(n) * diag(A);
}

func det(A)
{
   /* DOCUMENT det(A)
        Determinant of A.
    */

   if (dimsof(A)(2)==2 && dimsof(A)(3)==2) {
      return A(1,1)*A(2,2) - A(1,2)*A(2,1);
   } else {
      error, "not yet implemented";
   }
}

func spur(A)
{
   /* DOCUMENT spur(A)
        Returns the trace of matrix A, tr A.
    */
   for (i=1,s=0; i<=dimsof(A)(2); i++) {
      s += A(i,i);
   }
   return s;
}

func matrix_cut(&A, beg, end)
{
   /* DOCUMENT matrix_cut(A, beg, end)
        Reshapes the matrix A in situ to A(1+beg:N-end,1+beg:N-end).
    */

   if (is_void(beg) || is_void(end) || is_void(A)) error, "not enough parameters";

   N = dimsof(A)(3);

   if (dimsof(A)(1)==2) A = A(1+beg:N-end,1+beg:N-end);
   else if (dimsof(A)(1)==3) {
      B = array(double, dimsof(A)(2), dimsof(A)(3)-beg-end, dimsof(A)(4)-beg-end);
      for (i=1; i<=dimsof(A)(2); i++)  B(i,,) = A(i,1+beg:N-end,1+beg:N-end);
      A = B;
   }

}

func corr_coeff(A, zero=)
{
   /* DOCUMENT corr_coeff(A, zero=)
        Returns the matrix r with r(i,j) = A(i,j)/sqrt(A(i,i)*A(j,j)).
	If zero is given, sets any zero or negative diagonal element to this number.
    */
 
   if (!is_void(zero)) {
      w = where(diag(A)<=0);
      if (is_array(w))
	for (i=1; i<=dimsof(w)(2); i++)
	  A(w(i),w(i)) = zero;
   }
   return A/sqrt(diag(A)()(,-)*diag(A)(-,)); 
}

func ffeoc_rroc(r, sigma)
{
   /* DOCUMENT ffeoc_rroc(r, sigma)
        Inverse of corr_coeff.
    */

   N = dimsof(r)(2);
   if (N!=dimsof(r)(3) || N!=dimsof(sigma)(2)) error, "wrong dimensions";
   A = array(double, N, N);
   
   for (i=1; i<=N; i++) {
      for (j=1; j<=N; j++) {
         A(i,j) = sigma(i)*sigma(j)*r(i,j);
      }
   }

   return A;
}


func test_unity(a, b)
{
   /* DOCUMENT test_unity(a, b)
        Checks whether a*b = 1 (unit matrix). Sets array unity.
        If b not given, checks whether a is unity.
    */

   extern unity;

   if (!is_void(b)) {
      unity = a(,+)*b(+,);
   } else unity = a;

   d = diag(unity);

   write, format="diag: % .4e ... % .4e     ", min(d), max(d);
   for (u=unity,i=1; i<=dimsof(unity)(2); i++) u(i,i) = 0.0;
   write, format="off-diag: % .4e ... % .4e     ", min(u), max(u);

   write, format="cond = %.2e\n", LUrcond(a);
}

func test_symmetry(a, niceoutput=)
{
   /* DOCUMENT test_symmetry(a, niceoutput=)
         Checks whether a is symmetric.
    */

   maxd = -1.e30;
   maxq = -1.e30;
   eps  = 1.e-50;
   ljq = ljd = [0,0];

   n = dimsof(a)(2);

   for (i=1; i<=n; i++) {
      for (j=i; j<=n; j++) {
	 d = a(i,j) - a(j,i);
	 if (d>maxd) {
	    maxd = d;
	    ijd  = [i,j];
	 }
	 if (abs(a(j,i))>eps) {
            q = a(i,j)/a(j,i);
            if (q>maxq) {
               maxq = q;
               ijq  = [i,j];
            }
         }
      }
   }
   if (niceoutput==0) {
      write, format="max diff=%e (%d, %d), max quot=%e (%d %d)\n", maxd,
         ijd(1), ijd(2, maxq, ijq(1), ijq(2));
    } else {
        thresh = 1.0e-10;
        if (abs(maxd)<thresh && abs(maxq-1.0)<thresh)
	    write, format="%s", "matrix is symmetric\n";
        else
              write, format="%s", "matrix is not symmetric\n";
    }
}

func write_vector(v, file=, format=, tex=)
{
   /* DOCUMENT write_matrix(A, file=, format=, tex=)
        Writes matrix A  in rows format.
        If file not given, writes to stdout.
      SEE ALSO: write_matrix, write_matrix_col
    */

   if (is_void(tex)) tex = 0;

   N = dimsof(v)(2);
   if (is_void(format)) format = "% .6e";
   form = format + " ";

   for (i=1; i<=N; i++) {
      write, file, format=form, double(v(i));
      if (tex && i<N) write, file, format="%s", "& ";
      if (tex && i==N) write, file, format="%s", "\\\\";
   }
   writeNL, file;
}


func write_matrix(A, file=, format=, tex=)
{
   /* DOCUMENT write_matrix(A, file=, format=, tex=)
        Writes matrix A  in rows format.
        If file not given, writes to stdout.
      SEE ALSO: write_matrix_col, write_vector
    */

   if (is_void(tex)) tex = 0;
   
   N1 = dimsof(A)(2);
   N2 = dimsof(A)(3);
   if (is_void(format)) format = "% .6e";
   form = format + " ";
   
   for (i=1; i<=N1; i++) {
      for (j=1; j<=N2; j++) {
         write, file, format=form, double(A(i,j));
         if (tex && j<N2) write, file, format="%s", "& ";
         if (tex && i<N1 && j==N2) write, file, format="%s", "\\\\";
      }
      writeNL, file;
   }
}

func write_matrix_col(A, x, file=, format=, formatx=)
{
   /* DOCUMENT write_matrix_col(A, x, file=, format=)
        Writes square matrix A in column format with index array x.
        If file not given, writes to stdout.
      SEE ALSO: write_matrix, write_vector
    */

   N = dimsof(A)(2);
   if (dimsof(x)(2)!=N || dimsof(A)(3)!=N) error, "wrong dimensions";

   if (is_void(format)) format = "% .6e";
   if (is_void(formatx)) formatx = format;
   form = formatx + " " + formatx + " " + format + "\n";
   
   x1 = x()(-:1:N,);
   x2 = x()(,-:1:N);
   write, file, format=form, x1, x2, A;
}

func read_matrix(name, head=)
{
   /* DOCUMENT read_matrix(name, head=)
        Returns the matrix stored in file "name", e.g. created by write_matrix.
        head indicates the number of header lines. NEWMATORD
    */

   write, format="%s\n", "read_matrix NEWMATORD";

   if (is_void(head)) head = 0 ;

   N = numberoflines(name);
   A = array(double, N, N);

   f = open(name);
   for (i=1; i<=head; i++) rdline, f;   /* header */
   read, f, A;
   close, f;

   return transpose(A);
}


/* ============================================================ *
 * Statistics.                                                  *
 * ============================================================ */

func variance(x)
{
   /* DOCUMENT rms(x)
         Returns the (optimal, biased estimator of the) rms of x,
         where x can be an arbitrary-dimensional array.
	 The unbiased estimator is rms_unb = sqrt(n/(n-1))*variance(x)
	 where n = dimsof(x(*)).
    */

   return x(*)(rms);
}

func test_variance(x)
{
   n = dimsof(x)(2);

   extern a, s;

   for (a=0.0,i=1; i<=n; i++) {
      for (j=1; j<=n; j++) {
	 a += x(i,j);
      }
   }
   a = a/n^2;

   for (s=0,i=1; i<=n; i++) {
      for (j=1; j<=n; j++) {
	 s += (a-x(i,j))^2;
      }
   }
   sopt = sqrt(s/n^2);
   sunb = sqrt(s/(n^2-1.));

   write, format="avg, rms_opt, rms_unb=%e, %e, %e\n", a, sopt, sunb;
}

func bootstrapping(nstart, nboot)
{
   /* DOCUMENT bootstrapping(nstart, nboot)
        Draws nboot pairwise different entries out of [1:nstart] and
        returns them as an array.
    */

   if (nstart<nboot) error("nstart<nboot");

   istart = indgen(nstart);
   iend   = array(int, nboot);

   for (i=1; i<=nboot; i++) {

      nact = nstart-i+1;
      bt = int(random()*nact+1);
      bt = istart(bt);
      iend(i) = bt;

      if (i==nboot) break;

      ibt = where(istart==bt)(1);
      if (ibt==1) istart = istart(2:nact);
      else if (ibt==nact) istart = istart(1:nact-1);
      else istart = grow(istart(1:ibt-1), istart(ibt+1:nact));
 
   }

   return iend;
}

func real_bootstrapping(n)
{
   /* DOCUMENT real_bootstrapping(n)
        Returns n random numbers between 1 and n. Duplicates possible.
        This is the `real' bootstrapping.
    */

   return int(random(n)*n+1);
}

func random_dis(Nmax, N)
{
   /* DOCUMENT random_dis(N, Nmax)
         Like bootstrapping, but faster for large numbers. Returns sorted list.
    */

   extern xs;
   x = int(random(N)*Nmax + 1);
   xs = x(sort(x));

   do {

      for (i=1,d=0; i<=N-1; i++) {
         if (xs(i)==xs(i+1)) {
            d++;
            if (i==N-1) xs= grow(xs(1:N-1), int(random()*N + 1));
            else xs = grow(xs(1:i), xs(i+2:N), int(random()*N + 1));
         }
      }
      write, format="%d doubles\n", d;
      xs = xs(sort(xs));
   } while (d>0);

   return xs;
}

func random_init(noreread=)
{
   /* DOCUMENT random_init(void)
        Initializes the random seed by the current time (seconds).
        If noreread==1, reads the file .randomseed (must exist) for the seed.
      NOTE: Obsolete! And does not work (on Mac?).
      SEE ALSO: randomize
    */

   seed = array(int);
   if (is_void(noreread) || noreread==0) {
      $date +%N > .randomseed;
   }
   f = open(".randomseed");
   read, f, seed;
   close, f;
   //$rm .randomseed;

   seed = seed/1e9;
   /* write, format="initialize random generator with %f\n", seed; */
   random_seed, seed;
}

/* Needed by gasdev */
iset=0;
extern gset;
 
func gasdev(void)
{
   /* DOCUMENT gasdev
        gasdev returns a random variable of a Gaussian distribution N(0, 1).
        See Numerical Recipes.
    */

   a = gset;           /* for some obscure reason, a statement invoking gest is needed */
   if (iset==0) {
      do {
         v1 = 2.0*random() - 1.0;
         v2 = 2.0*random() - 1.0;
         rsq = v1^2 + v2^2;
      } while(rsq >= 1.0 || rsq==0.0);
      fac = sqrt(-2.0*log(rsq)/rsq);
      gset = v1*fac;
      iset = 1;
      return v2*fac;
   } else {
      iset = 0;
      return gset;
   }
}

func gastest(void)
{
   /* DOCUMENT gastest(void)
        Test for gasdev. The height of the resulting Gaussian histogram
        depends on N and on the binsize.
    */

   extern x;
   N = 50000;

   x = array(double, N);
   for (i=1; i<=N; i++) {
      x(i) = gasdev();
   }

   f = open("hist", "w");
   write, f, format="%.7f\n", x;
   close, f;
}

func sobseq(void)
{

}


/* ============================================================ *
 * FFT related functions.                                       *
 * ============================================================ */

func prefactor(N, Delta)
{
   /* DOCUMENT prefactor(N, Delta)
        Conversion between pixels and frequencies for fft. N is the
        length of the array, Delta the spacing between two pixels, as
        interpreted physically. For some index i, the physical
        quantity (frequency) is u_i = i*prefactor(N, Delta).
    */

   return twopi/(Delta*N);
}

func i2Ff(i, N)
{
   /* DOCUMENT i2Ff(i, N)
        i2Ff = index to Fourier frequency
        Returns the Fourier frequency f according to the index i.
        Can also be used in real space, for an array which is to be fft'd.
        f has to be multiplied by the proper spacing (Delta in real space,
        prefactor = 2*pi/(N*Delta) in Fourier space) to get the
        physical frequency.
    */

   if (i<=N/2) return i-1;
   return i-N-1;
}

func frequencies(n, Delta)
{
   /* DOCUMENT frequencies(N, Delta)
         Creates an array of length N with Fourier frequencies
         corresponding to a spacing in real space of Delta. The
         Fourier spacing is 1/(N*Delta).
    */

   fp = span(0, 1./(2*Delta), n/2+1);
   fm = -fp(n/2:2:-1);
   ff = array(float, n);
   ff(1:n/2+1) = fp;
   ff(n/2+2:n) = fm;

   return ff
}

func iscale(i, N, Ny, F)
{
   /* DOCUMENT iscale(i, N, Ny, F)
        Returns the index if an array of length N is to be expanded to length Ny.
        Used by expandField.
    */

   if (F==1) {
      if (i<=N/2) return i;
      return Ny-N+i;
   } else if (F==0) {
      return i + (Ny-N)/2;
   } else error("wrong F");
}

func expandField(x, N, expand, F)
{
   /* DOCUMENT expandField(x, N, expand)
        Expands the real 2D-array x (NxN) into a bigger array y
        (N*expand, N*expand). Returns y.
	F=0: x is centered in new array.
	F=1: x is split to the edges of the new array (old method).
    */

   if (is_void(F)) F = 0;

   Ny = N*expand;

   y = array(float, Ny, Ny);

   for (i=1; i<=N; i++) {
      iy = iscale(i, N, Ny, F);

      for (j=1; j<=N; j++) {
         jy = iscale(j, N, Ny, F);
         y(iy,jy) = x(i,j);
      }

   }

   return y;
}

func smooth1d(a, n, factor, scale=)
{
   /* DOCUMENT smooth1d(a, n, factor, scale=)
        Smoothes the 1d-array a with a (Gaussian) kernel of width
	n/factor.
	If scale>1 rescales the array to reduce aliasing.
    */

   if (is_void(scale)) scale = 1;

   N = n*scale;
   gauss = array(0.0, N);
   Sigma = [n/double(factor)];

   for (i=1; i<=N; i++) {
      gauss(i) = exp(mvdens_log_pdf([i], [0], Sigma));
   }

   ascale = array(0.0, N);
   ascale(n*(scale-1)/2+1:n*(scale-1)/2+n) = a;
   
   A     = fft(ascale);
   Gauss = fft(gauss).re;
   Res   = A*Gauss;
   res   = fft(Res, -1)/N;

   res   = res(n*(scale-1)/2+1:n*(scale-1)/2+n);
   
   return res.re;
}

func smooth2d(a, n, factor=, Sigma=, scale=)
{
   /* DOCUMENT smooth2d(a, n, factor=, Sigma=, scale=)
        Smoothes the 2d-array a with a (Gaussian) kernel of width
	[n(1)/factor, n(2)/factor];
	If scale>1 rescales the array to reduce aliasing.
    */

   if (is_scalar(n)) n = [n,n];
   if (is_void(scale)) scale = 1;

   if (is_void(factor) + is_void(Sigma) != 1) error, "One and only one of factor or Sigma have to be specified";

   N = n*scale;
   gauss = array(0.0, N(1), N(2));
   if (!is_void(factor)) {
      Sigma = array(0.0, 2, 2);
      Sigma(1,1) = n(1)/double(factor);
      Sigma(2,2) = n(2)/double(factor);
   }

   for (i=1; i<=N(1); i++) {
      for (j=1; j<=N(2); j++) {
	 x = [i,j];
	 gauss(i,j) = exp(mvdens_log_pdf(x, [0,0], Sigma));
      }
   }

   ascale = array(0.0, N(1), N(2));
   ascale(n(1)*(scale-1)/2+1:n(1)*(scale-1)/2+n(1), n(2)*(scale-1)/2+1:n(2)*(scale-1)/2+n(2)) = a;
   
   A     = fft(ascale);
   Gauss = fft(gauss).re;
   Res   = A*Gauss;
   res   = fft(Res, -1)/N(1)/N(2);

   res   = res(n(1)*(scale-1)/2+1:n(1)*(scale-1)/2+n(1), n(2)*(scale-1)/2+1:n(2)*(scale-1)/2+n(2));
   
   return res.re;
}


/* ============================================================ *
 * Vector manipulation.                                         *
 * ============================================================ */

func sortindex(x)
{
   /* DOCUMENT sortindex(x)
         Returns the array x, sorted in ascending order
    */

   xs = sort(x);
   return x(xs);
}

func perm(p)
{
   /* DOCUMENT perm(p)
        Returns the permutation #p of (1,2,3).
        For p<0 returns odd permutations.
    */

   if (p>=0) return ([0,1,2] + p%3)%3 + 1;
   else return ([0,2,1] - p%3)%3 + 1;
}

func normalize(x)
{
   /* DOCUMENT normalize(x)
         Returns x/||x||.
    */

   return x/sqrt(sum(x^2));
}

/* ============================================================ *
 * Maths Misc.                                                  *
 * ============================================================ */

func round_n(x, n=)
{
   /* DOCUMENT round_n(x, n=)
        Rounds the floating point variable x (to n digits after the
        point if n given.
    */

   if (is_void(n)) n = 0;

   return sign(x)*(int(abs(x)*10^n+0.5))/10.0^n;
}

func interp2d(f, p, delta)
{
   /* DOCUMENT interp2d(f, p, delta)
        Returns f interpolated for position p+delta where p is the 2d-index and
        delta \in [0,1[^2.
    */

   return delta(1)*delta(2)*f(p(1),p(2)) \
     + delta(1)*(1.0-delta(2))*f(p(1),p(2)+1) \
     + (1.0-delta(1))*delta(2)*f(p(1)+1,p(2)) \
     + (1.0-delta(1))*(1.0-delta(2))*f(p(1)+1,p(2)+1);
}

func prod(x)
{
   /* DOCUMENT prod(x)
        Return the product of all elements of the vector x,
        Prod_i=1^n x(i)
    */

   p = 1;
   n = dimsof(x)(2);
   for (i=1; i<=n; i++) {
      p *= x(i);
   }

   return p;
}

/* ============================================================ *
 * Eigensystems.                                                *
 * ============================================================ */

func ROTATE(a,i,j,k,l,s,tau)
{
   /* DOCUMENT ROTATE(a,i,j,k,l,s,tau)
        Used by jacobi.
      SEE ALSO: jaobi
    */

   g = a(i,j);
   h = a(k,l);
   a(i,j) = g - s*(h+g*tau);
   a(k,l) = h + s*(g-h*tau);
}

func sum_upper_right(a)
{
   /* DOCUMENT sum_upper_right(a)
        Returns the sum of the elements above the diagonal. Used by jacobi.
    */

   n = dimsof(a)(2);

   for (s=0,i=1; i<=n-1; i++)
     for (j=i+1; j<=n; j++)
       s += a(i,j);

   return s;
}

func choldc(&a, &p)
{
   /* DOCUMENT choldc(&a, &p)
        Cholesky decomposition of a, NR p.97. Returns 1 if the
        algorithm failes which means that a is not positive
        definite. Returns 0 upon success.
   */

   n = dimsof(*a)(2);

   for (i=1; i<=n; i++) {
      for (j=i; j<=n; j++) {
         for (s=(*a)(i,j),k=i-1; k>=1; k--) {
            s -= (*a)(i,k)*(*a)(j,k);
         }
         if (i==j) {
            if (s<=0.0) {
               /* matrix not positive definite */
               return 1;
            }

            (*p)(i) = sqrt(s);
         } else (*a)(j,i) = s/(*p)(i);
      }
   }

   /* decomposition successful */
   return 0;
}

func cholesky(a, &l)
{
   /* DOCUMENT cholesky(a)
        Calculates l where l l^t = a. Returns 0 upon success and 1 if
        a is not positive definite.
    */

   n = dimsof(a)(2);
   p = array(double, n);
   L = array(double, n, n);
   A = a;
   if (choldc(&A,&p)==1) return 1;

   for (i=1; i<=n; i++) {
      for (j=1; j<i; j++) {
         L(i,j) = A(i,j);
      }
      L(i,i) = p(i);
   }

   l = L;
   return 0;
}

func jacobi(input, iter)
{
   /* DOCUMENT jacobi(input,iter)
        Calculates the eigenvalues and -vector of the symmetric matrix input.
        Sets the extern variables eigenvalue and eigenvector.
	eigenvector(i,j) = i-th component of j-th eigenvector (so the
        eigenvectors are the columns of this matrix).
        With NEWMATORD, the consistency check can be done as
        input(,+)*eigenvector(+,i)/eigenvector(+,i) != eigenvalue(i).
        Default: iter=50
      SEE ALSO: eigenout, ROTATE
    */

   if (is_void(iter)) iter = 50;

   a = input;
   n = dimsof(a)(2);
   v = unit(n,n);

   z = array(0.0, n);
   b = d = diag(a);
   
   nrot=0;
   for (i=1; i<=iter; i++) {

      if (DEBUG) write, format="%c", '.';

      sm = sum_upper_right(a);

      if (sm==0.0) goto end_jacobi;

      if (i<4) tresh = 0.2*sm/n^2;
      else tresh = 0.0;

      for (ip=1; ip<=n-1; ip++) {
         for (iq=ip+1; iq<=n; iq++) {
            g = 100.0*abs(a(ip,iq));

            if (i>4 && abs(d(ip))+g==abs(d(ip))
                && abs(d(iq))+g==abs(d(iq))) {
               a(ip,iq) = 0.0;
            } else if (abs(a(ip,iq))>tresh) {
               h = d(iq) - d(ip);

               if (abs(h)+g==abs(h))
                 t = a(ip,iq)/h;
               else {
                  theta = 0.5*h/a(ip,iq);
                  t = 1.0/(abs(theta)+sqrt(1.0+theta^2));
                  if (theta<0.0) t = -t;
               }

               c = 1.0/sqrt(1+t^2);
               s = t*c;
               tau = s/(1.0+c);
               h = t*a(ip,iq);
               z(ip) -= h;
               z(iq) += h;
               d(ip) -= h;
               d(iq) += h;
               a(ip,iq) = 0.0;

               for (j=1; j<=ip-1; j++) nn = ROTATE(a,j,ip,j,iq,s,tau);
               for (j=ip+1; j<=iq-1; j++) nn = ROTATE(a,ip,j,j,iq,s,tau);
               for (j=iq+1; j<=n; j++) nn = ROTATE(a,ip,j,iq,j,s,tau);
               for (j=1; j<=n; j++) nn = ROTATE(v,j,ip,j,iq,s,tau);
               ++nrot;
            }
         }
      }

      b = b + z;
      d = b;
      z = array(0.0, n);
   }
   error("Too many iterations in routine jacobi");

 end_jacobi:

   if (DEBUG) writeNL;

   extern eigenvalue, eigenvector;

   eigenvalue  = d;
   eigenvector = v;
}

func eigenout(file)
{
   /* DOCUMENT eigenout(file)
        Prints the (sorted) eigenvalues and eigenvectors, as in ~/bin/fisher-variance.
	If file not given, writes to stdout.
    */

   n = dimsof(eigenvalue)(2);
   s = sort(eigenvalue);
   //s = indgen(n);

   if (is_array(where(eigenvalue<0))) write, file, format="%s\n", "WARNING: eigenvalue<0"; 
   write, file, format="%s", "ew="
   for (area=pi,i=1; i<=n; i++) {
      write, file, format="%.3e ", eigenvalue(s(i));
      area *= sqrt(abs(eigenvalue(i)));
   }
   write, file, format=" Area=%.5f\n\n", area;

   for (i=1; i<=n; i++) {
      for (j=1; j<=n; j++)
        write, file, format="%9.3f ", eigenvector(i,s(j));
      write, file, "";                                       /* newline */
   }
}

/* ============================================================ *
 * Plotting.                                                    *
 * ============================================================ */

extern plot_type;
extern plot_color;
plot_type = 1;
plot_color = 1;

func plot_reset(l=)
{
   if (is_void(l)) l = 1;

   extern plot_type;
   plot_type = l;

   extern plot_color;
   plot_color = 1;
}

func plot(y, x, black=, legend=, hide=, type=, width=, color=, closed=, smooth=, marks=, marker=, mspace=,
          mphase=, rays=, arrowl=, arroww=, rspace=, rphase=)
{
   /* DOCUMENT
        Calls plg without mark (default) and with successive line types (default).
    */

   if (is_void(black)) black = 1;
   
   if (!is_void(marks)) m = marks;
   else m = 0;

   extern plot_type;
   if (!is_void(type)) l = type;
   else l = plot_type;

   extern plot_color;
   if (!is_void(color)) c = color;
   else c = Color(,plot_color%dimsof(Color)(3));

   plg, y, x, legend=legend, hide=hide, type=l, width=width, color=c, closed=closed,
     smooth=smooth, marks=m, marker=marker, mspace=mspace, mphase=mphase, rays=rays,
     arrowl=arrowl, arroww=arroww, rspace=rspace, rphase=rphase;

   plot_type += 1;
   plot_color += 1;
}

func plgw(y, x, marks=, type=, width=, color=)
{
   /* DOCUMENT plgw(y, x, marks=, type=, width=, color=)
        Plots y(where(y)) against x(where(y))
    */

   plg, y(where(y)), x(where(y)), marks=marks, type=type, width=width,
     color=color;
}

func plot_errorbars(x, y, err, onlyerr=, type=, color=, width=, symbol=)
{
   /* DOCUMENT plot_errorbars(x, y, err, onlyerr=, type=, color=, width=, symbol=)
        Plots the function y versus x with errorbars err.
        If onlyerr=1, don't connect the points with lines.
      SEE ALSO: plp (yutils/plot.i)
    */

   if (is_void(onlyerr)) onlyerr = 0;
   
   if (onlyerr==0) plg, y, x, marks=0, type=type, color=color;
   if (!is_void(symbol)) plmk, y, x, marker=symbol, color=color;
   pldj, x, y-err/2, x, y+err/2, color=color, width=width;
}

func plot_line(a, b)
{
   /* DOCUMENT plot_line(a, b)
        Plots a line from a to y (both 2-d vectors).
    */

   plg, [a(2),b(2)], [a(1),b(1)], marks=0;
}

func plotcontour(Cov2, Theta, r, ownlevs, levdata, marks=, color=)
{
   /* DOCUMENT plotcontour(Cov2, Theta, r, ownlevs, levdata)
        Makes a contour plot of the 2-d array Cov2. Theta defines the
        grid. If ownlevs=1, the levels are determined automatically
        (by plc). If ownlevs=0, levdata = [lmin, lmax, Nl] is taken to
        create levels. If r=1, the correlation coefficient
        r_ij=Cov2_ij/sqrt(Cov2_ii Cov2_jj) is plotted.
        Defaults: Cov2=cov2, Theta=theta, r=0, ownlevs=0, levdata =
        [1e-13, 1e-10, 20].
        
    */

   if (is_void(Cov2)) Cov2 = cov2;
   if (is_void(Theta)) {
      if (is_void(theta)) {
	 Theta = indgen(dimsof(Cov2)(2));
      } else {
	 Theta = theta;
      }
   }
   if (is_void(r)) r = 0;
   if (is_void(ownlevs)) ownlevs = 0;
   if (is_void(levdata)) levdata = [1e-13, 1e-10, 20];
   if (is_void(marks)) marks=0;
   if (is_void(color)) color = Black;

   Ntheta = dimsof(Theta)(2);
   gx = Theta()(-:1:Ntheta,);
   gy = Theta()(,-:1:Ntheta);

   if (r==0) {

      levels    = spanl(levdata(1), levdata(2), int(levdata(3)));
      levelsneg = spanl(-levdata(1),-levdata(2), int(levdata(3)));

   } else {

      /* check !!! ??? !!! */
      levels    = span(0, 1, 20);
      levelsneg = span(-1, 0, 20);
      Cov2      = corr_coeff(Cov2);

   }

   if (ownlevs==0) {
      plc, Cov2, gy, gx, levs=levels, marks=marks, color=color;
      plc, Cov2, gy, gx, levs=levelsneg, marks=marks, type="dash", color=color;
   } else {
      plc, Cov2, gy, gx, marks=marks, color=color;
   }
}

func plot_corrcoeff(A, x, j, linetype=)
{
   /* DOCUMENT plot_corrcoff(A, x, j, type=linetype)
        Plots the correlation coefficients of the matrix A.
    */

   if (is_void(j)) j = 0;

   n = dimsof(A)(2);
   if (n!=dimsof(x)(2)) error("wrong dimensions");
   A = corr_coeff(A);

   for (i=1; i<=n; i++) {
      if (j>0 && i!=j) continue;
      plg, A(i,), x, marks=0, type=linetype;
      str = swrite(format="#%d", i);
      if (j==-1) read, prompt=str, format="\n", xx;
   }

}

func ellipse(mean, var, sigma=, N=)
{
   /* DOCUMENT ellipse(mean, var, sigma=, N=)
        
    */

   if (is_void(sigma)) sigma = 1.0;
   
   if (is_void(N)) N = 100;

   jacobi, var;
   if (eigenvalue(1)<=0 || eigenvalue(2)<=0) {
      //write, format="Warning: eigenvalue (%g,%g) < 0\n", eigenvalue(1), eigenvalue(2);
      return [];
   }
   ewsqrt = sqrt(eigenvalue)*sigma;

   X = array(double, 2, N);
   Z = array(double, 2, 2*N);
   X(1,) = span(-ewsqrt(2), ewsqrt(2), N);
   X(2,) = sqrt(1-(X(1,)/ewsqrt(2))^2)*ewsqrt(1);

   cospsi = eigenvector(1,2);
   sinpsi = eigenvector(2,2);
   Z(1,1:N) = cospsi*X(1,) - sinpsi*X(2,);
   Z(2,1:N) = sinpsi*X(1,) + cospsi*X(2,);

   Z(1,N+1:2*N) = cospsi*X(1,N:1:-1) + sinpsi*X(2,N:1:-1);
   Z(2,N+1:2*N) = sinpsi*X(1,N:1:-1) - cospsi*X(2,N:1:-1);
   Z = Z + mean;

   return Z;
}

func draw_ellipse(Z, delta=, type=, width=, color=)
{
   /* 1-Sigma-Ellipse. MKDEBUG: This function does not *
    * what it's supposed to do. delta is not used.     */

   if (is_void(delta)) delta = [0.05, 0.05];
   if (is_void(type)) type = "solid";
   if (is_void(width)) width = 1;
   if (is_void(color)) color = Black;

   if (!is_void(Z)) {
      plg, Z(2,), Z(1,), marks=0, type=type, color=color, width=width;
   }
}

func colorbar(cmin, cmax)
{
   /* DOCUMENT colorbar
               colorbar, cmin, cmax
        draw a color bar to the right of the plot.  If CMIN and CMAX
        are specified, label the top and bottom of the bar with those
        numbers.
    */

  plsys, 0;
  pli, span(0,1,200)(-,), .625,.46,.67,.84, legend="";
  plg, [.46,.84,.84,.46],[.67,.67,.625,.625], closed=1,
    marks=0,color="fg",width=1,type=1,legend="";

  plsys, 1;  /* assumes there is only one coordinate system */

  if (!is_void(cmin)) {
    plt, pr1(cmin), .6475,.46, justify="CT";
    plt, pr1(cmax), .6475,.84, justify="CB";
  }
}


func colorbartop(c, pmin=, pmax=, height=)
{
   /* DOCUMENT colorbartop
               colorbartop, c, pmin=, pmax=, height=
        Draws a color on top of the plot. For values of the array c,
        labels are drawn. The color range is the full range of the
        palette, or between [pmin; pmax]*100 percent of the range.
        c can be array of string or number.
        Default: pmin=0.0, pmax=1.0.

        Example: image for values between levmin and levmax is
        plotted.
        Logarithmic levels:
           c1 = 1e-5; c2 = 1e4; nc = 10;
           pmin = (log10(c1)-log10(levmin))/(log10(levmax)-log10(levmin));
           pmax = (log10(c2)-log10(levmin))/(log10(levmax)-log10(levmin));
           colorbartop, spanl(c1, c2, nc), pmin=pmin, pmax=pmax;

        Linear levels:
           c1 = 0; c2 = 9; nc = 10;
           pmin = (c1-levmin)/(levmax-levmin);
           pmax = (c2-levmin)/(levmax-levmin);
           colorbartop, span(c1, c2, nc), pmin=pmin, pmax=pmax;
    */

  if (is_void(height)) height = 18;
  port = viewport();
  /* dx   = [0.03, 0.01]; */
  dx   = [0.0, 0.01];
  width = .035;
  xll  = port(1:4:3)+dx;                        /* lower left corner */
  xur  = [port(2)-dx(1), port(4)+dx(2)+width];
  /* upper right corner of colorbar */

  plsys, 0;

  if (is_void(pmin)) pmin = 0.0;
  if (is_void(pmax)) pmax = 1.0;
  pli, bytscl(span(pmin,pmax,100), cmin=0, cmax=1)(,-), xll(1), xll(2),
    xur(1), xur(2), legend="";
  plg, [xll(2), xur(2), xur(2), xll(2)], [xll(1), xll(1), xur(1), xur(1)],
    closed=1,marks=0,color="fg",width=1,type=1,legend="";

  plsys, 1;  /* assumes there is only one coordinate system */

  if (!is_void(c)) {
     if (is_scalar(c)) error, "colorbar labels not given as array";
     nc = dimsof(c)(2);

     x = xur(1) - xll(1);
     for (i=1; i<=nc; i++) {
        xi = double(i-1)/(nc-1)*x + xll(1);
        y  = xur(2) + 0.005;

        if (wtype(c(i))==T_STRING) text = c(i);
        else text = pr1(c(i));

        if (i==1) plt, text, xi, y, justify="LB", height=height;
        else if (i==nc) plt, text, xi, y, justify="RB", height=height;
        else plt, text, xi, y, justify="CB", height=height;
     }
  }
}

func pltitle2(title, deltay)
/* DOCUMENT pltitle2, title, deltay=
     Plot TITLE centered above the coordinate system for any of the
     standard Gist styles.  You may want to customize this for other
     plot styles.
   SEE ALSO: plt, xytitles, pltitle
 */
{
   if (is_void(deltay)) deltay = 0.0;
   
  port= viewport();
  plt, title, port(zcen:1:2)(1), port(4)+0.02+deltay,
    font=pltitle_font, justify="CB", height=pltitle_height;
}

func pltextspecial(t, p, height=, logx=, logy=, font1=, font2=, epsy=, justify=, tosys=)
{
   /* DOCUMENT pltextspecial(t,p,height=)
        Plots t(1)_t(2) at position [p(1),p(2)].
    */

   if (is_void(height)) height=14.0;
   if (is_void(logx)) logx=0;
   if (is_void(logy)) logy=0;
   if (is_void(epsy)) epsy=1;
   if (is_void(tosys)) tosys=0;

   l = limits();
   if (logx==0) Lx = l(2)-l(1);
   else         Lx = log(l(2))-log(l(1));
   if (logy==0) Ly = l(4)-l(3);
   else         Ly = log(l(4))-log(l(3));

   pt = [Lx,Ly]/450.;
   dp = [height^0.65*1.5*pt(1)*epsy, height/4.*pt(2)];

   p2 = array(double, 2);
   if (loxy==0) p2(1) = p(1) + dp(1);
   else         p2(1) = p(1)*exp(dp(1));
   if (logy==0) p2(2) = p(2) - dp(2);
   else         p2(2) = p(2)/exp(dp(2));

   plt, t(1), p(1), p(2), tosys=1, font=font1, justify=justify, tosys=tosys;
   if (t(2)!="") plt, t(2), p2(1), p2(2), tosys=1, height=height*0.75,
                   font=font2, justify=justify, tosys=tosys;
}

func plot_circle(pos, size=, sizey=, width=, color=)
{
   Ncirc = 100;
   if (is_void(size)) size = 1;
   if (is_void(sizey)) sizey = size;
   if (is_void(width)) width = 1;

   t = span(0, 2*pi, Ncirc);
   plg, sizey*cos(t)+pos(2), size*sin(t)+pos(1), marks=0, color=color, width=width;

}

func my_eps(file)
{
   eps, "tmp_yor.eps";
   str = swrite(format="add_comment_to_ps.pl -l y -o %s.eps tmp_yor.eps", file);
   system(str);
   system("rm tmp_yor.eps");
}


func my_is_scalar(x)
{
   /* DOCUMENT my_is_scalar(object)
         returns 1 if OBJECT is a scalar, else 0.
	 This function should be the same as std.i:is_scalar.
      SEE ALSO: is_scalar, is_array, is_func, is_void, is_range, is_struct, is_stream
    */
   return is_array(x) && !dimsof(x)(1);
}

extern T_CHAR, T_SHORT, T_INT, T_LONG, T_FLOAT, T_DOUBLE, T_STRING, T_POINTER;
T_CHAR = 0; T_SHORT = 1; T_INT = 2; T_LONG = 3; T_FLOAT = 4;
T_DOUBLE = 5; T_STRING = 6; T_POINTER = 7;

func wtype(topic)
{
   name = nameof(structof(topic));

   if (name=="char") typeID = T_CHAR;
   else if (name=="short") typeID = T_SHORT;
   else if (name=="int") typeID = T_INT;
   else if (name=="long") typeID = T_LONG;
   else if (name=="float") typeID = T_FLOAT;
   else if (name=="double") typeID = T_DOUBLE;
   else if (name=="compelx") typeID = T_COMPLEX;
   else if (name=="string") typeID = T_STRING;
   else if (name=="pointer") typeID = T_POINTER;
   else error, "unkown typeID: ", name;

   return typeID;
}

/* ============================================================ *
 * Input/Output.                                                *
 * ============================================================ */

func read_table(name, &table)
{
   /* DOCUMENT read_table(name, &table)
         Reads the file "name" and fills the ncol*nlines-dimensional array
	 table. Inspired by R.
    */

   nl    = numberoflines(name, nh);
   if (nl==-1) { table = void; return -1; }
   nc    = numberofcolumns(name, all=0);

   //write, format="nl, nc, nh = %d %d %d\n", nl, nc, nh;
   if (nc<=0 || nl<=0 || nl-nh<=0) {
      write, format="Error: could not read file (nl=%d, nc=%d, nh=%d)\n", nl, nc, nh;
      table = void;
      return -1;
   }

   table = array(double, nc, nl-nh);

   f = open(name);
   for (i=1; i<=nh; i++) { hd = rdline(f); }
   read, f, table;
   close, f;
   return 0;
}

func epsall(maxwin, dir=, pre=)
{
   /* DOCUMENT epsall(maxwin, dir=, pre=)
        Hardcopies the contents of window 0,1,...,maxwin into the files
        "dir/pre{0,1,...}".
        Defaults: maxwin=2, dir=".", pre="yorick".
    */

   if (is_void(maxwin)) maxwin = 2;
   if (is_void(dir)) dir = ".";
   if (is_void(pre)) pre = "yorick";

   for (i=0; i<=maxwin; i++) {
      window, i;
      name = swrite(format="%s/%s%d", dir, pre, i);
      eps, name;
   }
}

func read_tak(name, g)
{
   /* DOCUMENT read_tak(name)
         Reads one of Takashi Hamana's LCDM simulation "name".
         Returns kappa/gamma1/gamma2 (g=0/1/2) as a N x N - array (N=1024).
         Default: g=0
    */

   if (is_void(g)) g = 0;

   N = 1024;
   y = array(float, N, N);      /* data is single precision (f77 real)! */

   f = open(name, "rb");

   _read, f, 4+4*(N^2+2)*g, y;
   close, f;

   return y;
}

func read_tak_new(name, &res)
{
   f = open(name, "rb");
   ndat = array(int);

   //_read, f, 4, ndat;
   res = array(float, 8, 100);
   _read, f, 8, res;
   close, f;
}

func writeNL(file)
{
   write, file, format="%s", "\n";
}

func ls(dir, endung=)
{
   /* DOCUMENT ls(dir, endung=)
        Returns an array of strings with all files in directory dir.
        Default: dir="."
    */

   if (is_void(dir)) dir = ".";

   f = open(".yorick.ls.tmp", "w");
   if (is_void(endung)) write, f, format="ls -1 %s\n", dir;
   else write, f, format="ls -1 %s/*.%s\n", dir, endung;
   close, f;
   
   $source .yorick.ls.tmp > .yorick.ls.tmp2;

   n = numberoflines(".yorick.ls.tmp2");
   x = array(string, n);
   f = open(".yorick.ls.tmp2", "r");
   for (i=1; i<=n; i++) {
      x(i) = rdline(f);
   }
   close, f;
   $rm ".yorick.ls.tmp";
   $rm ".yorick.ls.tmp2";

   return x;
}

/* ============================================ *
 * Multi-variate Gaussian                       *
 * ============================================ */

/* mrvnorm obsolete, replace by mvdens */
struct mvrnorm {
  int ndim;
  double *pmean, *pvar;
}

func mvrnorm_read(name)
{
   /* DOCUMENT mvrnorm_read(name)
         Reads a file created by montecarlo.c:mvrnorm_print.
    */

   write, format="%s\n", "*** Warning: mvrnorm_read is obsolete ***";

   ndim  = numberoflines(name)-1;
   mean  = array(double, ndim);
   var   = array(double, ndim, ndim);

   f = open(name);
   read, f, mean;
   read, f, var;
   close, f;

   x = mvrnorm(ndim=ndim, pmean=&mean, pvar=&var);
   return x;
}

struct mvdens {
  int ndim, df;
  double *mean, *var;
}

struct mix_mvdens {
  int ncomp, ndim;
  mvdens *mvd;
  double *weight;
}

func mvdenspart_read(f)
{
   /* DOCUMENT mvdenspart_read(f)
        Reads a struct mvdens from file f.
    */

   header = array(int, 4);
   read, f, header;

   ndim = header(1);
   df   = header(2);
   mean = array(double, ndim);
   var  = array(double, ndim, ndim);

   read, f, mean;
   read, f, var;

   mvd = mvdens(ndim=ndim,df=df,mean=&mean,var=&var);
   return mvd;
}

func mvdens_read(name)
{
   /* DOCUMENT mvdens_read(name)
        Reads a file created by mvdens.c:mvdens_dump.
    */

   f = open(name);
   mvd = mvdenspart_read(f);
   close, f;
   return mvd;
}

func mix_mvdens_read(name)
{
   /* DOCUMENT mix_mvdens_read(name)
        Reads a file created by mvdens.c:mix_mvdens_dump.
    */

   f = open(name);
   header = array(int, 2);
   read, f, header;

   ncomp  = header(1);
   ndim   = header(2);
   mvd    = array(mvdens, ncomp);
   weight = array(double, ncomp);

   for (i=1; i<=ncomp; i++) {
      read, f, weight(i);
      mvd(i) = mvdenspart_read(f);
   }
   
   mix_mvd = mix_mvdens(ncomp=ncomp,ndim=ndim,mvd=&mvd,weight=&weight);
   return mix_mvd;
}

func mvdens_write(mvd, fname=)
{
   /* DOCUMENT mvdens_write(mvd, fname=)
        Writes a mvdens struct mvd to a file fname (default stdout).
    */

    if (is_void(fname)) f = stderr;
    else f = open(fname, "w");

    write, f, format="%d %d %d %d\n", mvd.ndim, mvd.df, mvd.ndim, 1;
    for (i=1; i<=mvd.ndim; i++) {
      write, f, format="%g ", (*(mvd.mean))(i);
    }
    writeNL, f;
    for (i=1; i<=mvd.ndim; i++) {
       for (j=1; j<=mvd.ndim; j++) {
          write, f, format="%g ",  (*(mvd.var))(i,j);
       }
       writeNL, f;		
   }
}

func mvdens_log_pdf(x, mu, Sigma)
{
   /* DOCUMENT mvGauss(x, mu, Sigma)
        Returns the log-density of a mvdens (at the moment only Gauss).
    */

   ndim = dimsof(x)(2);
   if (ndim!=dimsof(mu)(2) || ndim!=dimsof(Sigma)(2)) {
      error("Wrong dimensions");
      return;
   }

   d      = x-mu;
   if (dimsof(Sigma)(2)>1) {
      sigma  = LUsolve(Sigma);
      tmp    = sigma(,+)*d(+);
      e      = d(+)*tmp(+);
      /* TODO: Use cholesky decomp */
      logdet = log(det(Sigma));
   } else {
      e = d(1)^2/Sigma(1);
      logdet = log(Sigma(1));
   }

   r = -0.5*(ndim*log(2*pi) + e + logdet);

   return r;
}

func mvdens_draw(Sigma, quiet=)
{
   /* DOCUMENT mvdens_draw(Sigma, quiet=)
        Returns a random vector drawn from N(0, Sigma).
      SEE ALSO malmquist:rvGauss.
    */

   if (is_void(quiet)) quiet = 0;

   N = dimsof(Sigma)(2);
   u = array(0.0, N);

   for (i=1; i<=N; i++) {
      g = gasdev();
      u(i) = g;
   }

   i = cholesky(Sigma, L);
   if (quiet == 0) {
      write, stderr, format="Cholesky returns %d\n", i;
   }

   z = L(,+)*u(+);
   return z;
}



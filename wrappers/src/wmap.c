#include "wmap.h"

functions_wrapper_t *init_functions_CMB(error **err)
{
   functions_wrapper_t *init_func;

   //fprintf(stderr, "init_functions_CMB\n");

   init_func = init_functions_wrapper(read_from_config_CMB, init_CMB, likeli_CMB, NULL, deduced_CMB, print_CMB, err);
   forwardError(*err, __LINE__, NULL);

   return init_func;
}

void read_from_config_CMB(void **state, FILE *F, error **err)
{
   wmap_state *wmapstate;
   config_element c = {0, 0.0, ""};

   //fprintf(stderr, "read_from_config_CMB\n");

   wmapstate = malloc_err(sizeof(wmap_state), err);
   forwardError(*err, __LINE__,);

   CONFIG_READ_S(wmapstate, scamb_path, s, F, c, err);
   CONFIG_READ_S(wmapstate, data_path, s, F, c, err);
   CONFIG_READ_S(wmapstate, Cl_SZ_file, s, F, c, err);
   CONFIG_READ(wmapstate, lmax, i, F, c, err);
   CONFIG_READ(wmapstate, accurate, d, F, c, err);
   CONFIG_READ_S(wmapstate, model_file, s, F, c, err);

   *state = wmapstate;
}

int wmap_readline(int fds, char *buf)
{
  int i;

  i = 0;
  if (read(fds,buf,1)==0) {
    return -1;
  }

  while((buf[i]!='\n') && (buf[i]!='$')) {
    i++;
    if (read(fds,&(buf[i]),1)==0) {
      return -1;
    }
  }

  if (buf[i]=='$') {
    buf[i]='\0';
    return -1;
  }
  if (i==0) {
     return wmap_readline(fds,buf);
  }
  
  buf[i]='\0';
  //fprintf(stderr, "readline done, returning %d\n", i);

  return i;
}

void wmap_lkl_loop(void* data, error** er)
{
  char buf[2048];
  //double cls[(TT_MAX+1)*4];
  double *cls;
  double *clt,*cle,*clb,*clte;
  int err, lmax;
  double rr;
  forked_state *self;
  wmap_state *wmap;


#ifdef COMM_DEBUG
  fprintf(stderr, "wmap_lkl_loop\n");
#endif

  
  self = (forked_state*)data;
  wmap = (wmap_state*)(self->specific_data);

  //fprintf(stderr, "wmap_lkl_loop: Calling lkl_init with data_path='%s'\n", wmap->data_path);

  lkl_init(wmap->data_path, er);
  forwardError(*er, __LINE__,);

  sprintf(buf,"$");
  write(self->pout,buf,2);

  lmax = wmap->lmax;
  cls = malloc_err(sizeof(double)*(lmax+1)*4, er);
  forwardError(*er, __LINE__,);

  clt = cls;

  //cle = cls+TT_MAX+1;
  //clb=cls+(TT_MAX+1)*2;
  //clte=cls+(TT_MAX+1)*3;

  cle  = cls+lmax+1;
  clb  = cls+(lmax+1)*2;
  clte = cls+(lmax+1)*3;
  
  while(1) {
     err = read(self->pin, cls, sizeof(double)*(lmax+1)*4);
    if (err!=sizeof(double)*(lmax+1)*4) {
      rr = -1.0;
      sprintf(buf,"$\n$");
    } else {
      rr = compute_lkl(NULL, clt, clte, cle, clb, &lmax);

#ifdef COMM_DEBUG
      fprintf(stderr, "In wmap_lkl_loop: likelihood %g to be written to buffer\n", rr);
#endif

      sprintf(buf,"$%30.20g\n$",rr);
    }
    //fprintf(stderr, "writing ---%s---%g\n", buf, rr);
    write(self->pout,buf,strlen(buf));
  }

  free(cls);

#ifdef COMM_DEBUG
  fprintf(stderr, "Returning from wmap_lkl_loop\n");
#endif
}

int wmap_lkl_parent_init(void* data, error** err)
{
  forked_state* self;
  char buf[2048];

  self = (forked_state*)data;

#ifdef COMM_DEBUG
  int errV;
  fprintf(stderr, "wmap_lkl_parent_init: Calling 'sendAndCheck' for WMAP likelihood\n");
  errV =
#endif

  sendAndCheck(self,NULL,0,self->timout);

#ifdef COMM_DEBUG
  fprintf(stderr, "sendAndCheck WMAP likelihood, errV=%d (%s) (time-out %ld)\n", errV,
	  errV==0?"success":"failure", (long)self->timout);
#endif

  while (wmap_readline(self->pout,buf)>=0);

  pipe_cleanup(self->pout);

  return 1;
}

void init_CMB(common_like *like, error** err)
{
  wmap_state *wmapstate;
  char *tmp;
  FILE *FD;


#ifdef COMM_DEBUG
  fprintf(stderr, "init_CMB start\n");
#endif


  wmapstate = (wmap_state*)(like->state);

  if (wmapstate->lmax<2000) {
     fprintf(stderr, "Warning, lmax=%d is smaller than 2000, needed for 0.3%% precision on the Cl's (comparable with current CAMB precision)\n", wmapstate->lmax);
  }
  testErrorRetVA(wmapstate->lmax<1200, wmap_lmax, "lmax=%d is smaller then 1200 (WMAP7 likelihood lower limit)",
		 *err, __LINE__,, wmapstate->lmax);

#ifdef COMM_DEBUG
  struct timeval before, after, finally;
  gettimeofday(&before,NULL);
#endif
  
  /* Buffer for the angular power spectrum */
  wmapstate->cls = malloc(sizeof(double)*(wmapstate->lmax+1)*4);
  if (wmapstate->cls==NULL) {
    free(wmapstate);
    *err = addError(wmap_allocate,"Cannot allocate",*err,__LINE__);
    return;
  }
  
  /* Always set want_sigma_8 to 1 even if not needed */
  wmapstate->want_sigma_8 = 1;
  wmapstate->want_out_Cl  = 0;
  wmapstate->want_out_Pk  = 0;
  wmapstate->do_lensing   = 1;
  wmapstate->n_template   = 0;

  /* Default cosmological model */
  if (strcmp(wmapstate->model_file, "-")!=0) {
     FD = fopen_err(wmapstate->model_file, "r", err);
     forwardError(*err, __LINE__, );
     read_cosmological_parameters(&wmapstate->model, FD, err);
     forwardError(*err, __LINE__, );
     fclose(FD);
  } else {
     wmapstate->model = set_cosmological_parameters_to_default2(err);
     forwardError(*err, __LINE__, );
  }


  /* SZ power spectrum */
  if (strcmp(wmapstate->Cl_SZ_file, "-")==0) {
     wmapstate->template   = NULL;
  } else {
     wmapstate->n_template++;
     wmapstate->template = read_Cl_SZ_file(wmapstate->Cl_SZ_file, wmapstate->lmax, err);
     forwardError(*err, __LINE__, );
  }
  
  wmapstate->parout = malloc_err(sizeof(double)*(16+wmapstate->n_template), err);
  forwardError(*err, __LINE__, );

  /* Initializing camb */
#ifdef COMM_DEBUG  
  fprintf(stderr, "Initializing camb\n");
#endif

  tmp = wmapstate->scamb_path;
  wmapstate->camb = init_exec(1, wmapstate->scamb_path, &tmp, DEFAULT_TIMEOUT,
			 NULL, "camb", err);
  if (wmapstate->camb==NULL) {
    free(wmapstate);
    forwardError(*err,__LINE__,);
    // MKDEBUG: error not outside of bracket?
  }

#ifdef COMM_DEBUG
  fprintf(stderr, "Init camb: success\n");
#endif
  
#ifdef COMM_DEBUG
  gettimeofday(&after,NULL); 
#endif
  
  /* Initializing wmap likelihood */

#ifdef COMM_DEBUG
  fprintf(stderr, "Init wmap wmap_lkl_loop=%p\n", wmap_lkl_loop);
#endif

  wmapstate->lkl = init_forked(wmap_lkl_loop, wmapstate, DEFAULT_TIMEOUT, NULL, wmap_lkl_parent_init, "wmap_lkl", 0, err);
  if (wmapstate->lkl==NULL) {
    free(wmapstate);
  }
  forwardError(*err, __LINE__, );
  

#ifdef COMM_DEBUG  
  gettimeofday(&finally, NULL); 

  fprintf(stderr, "init_CMB: time spent: %d msec (camb %d, wmap_lkl %d)\n",
          abs((before.tv_sec-finally.tv_sec)*1000000+(before.tv_usec-finally.tv_usec))/1000,
          abs((after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec))/1000,
          abs((after.tv_sec-finally.tv_sec)*1000000+(after.tv_usec-finally.tv_usec))/1000);
#endif

#ifdef COMM_DEBUG
  fprintf(stderr, "init_CMB end\n");
#endif

}

/* ============================================================ *
 * Reads the SZ C_l's from a file, format:			*
 * l C_l							*
 * See Komatsu&Seljak (2002). For WMAP7, use the file           *
 * 'wmap_clsz_KKaQVW_l=2-10000_v4.txt' (V-band column) from     *
 * http://lambda.gsfc.nasa.gov.					*
 * ============================================================ */
double *read_Cl_SZ_file(const char *Cl_SZ_file, double lmax, error **err)
{
   double *cl;
   int ell;
   FILE *F;
   int l;

   F = fopen_err(Cl_SZ_file, "r", err);
   forwardError(*err, __LINE__, NULL);

   cl = malloc_err(sizeof(double)*(lmax+1), err);
   forwardError(*err, __LINE__, NULL);

   cl[0] = cl[1] = 0.0;
   l = 2;
   while (!feof(F) && l<=lmax) {
      fscanf(F, "%d %lg", &ell, cl+l);
      testErrorRetVA(l!=ell, wmap_ell, "Wrong l=%d found in file, expected %d", *err, __LINE__, NULL, ell, l);
      l++;
   }
   testErrorRetVA(ell<lmax, wmap_lmax, "Premature end of file at l=%d, you have to go to at least lmax+1=%g",
		  *err, __LINE__, NULL, ell, lmax+1);
   fclose(F);
   return cl;
}

void dealWithIssue(forked_state *frk, int err, int* retry, int* status, error **er)
{
  char buferr[256];
  int ret;
  
  ret = *retry;
  
  if ((err==DEAD) || (err==HANGUP)) {
    fprintf(stderr, "%s is dead or hung up (err=%d), trying again...\n", frk->name, err);
    *retry = 0;
    *status=WMAP_DEAD;
    launchMe(frk,er);
    forwardError(*er,__LINE__,);
    *status = WMAP_CHOKE;
    *retry = ret;
    if (ret == 0) {
      sprintf(buferr,"%s is dead (%d)",frk->name,err);
      *er = addError(wmap_choke,buferr,*er,__LINE__);  
    }
    return;
  }
  if ((err == WRITEERROR) || (err == READERROR)) {
    fprintf(stderr,"%s is in very bad shape (write or read error), aborting",frk->name);
    sprintf(buferr,"%s is in very bad shape (write or read error), aborting",frk->name);
    *status=WMAP_DEAD;
    *retry=0;
    *er = addError(wmap_dead,buferr,*er,__LINE__);
    return; 
  }
  if (err == TIMEOUT) {
    fprintf(stderr,"%s time outed",frk->name);
    sprintf(buferr,"%s time outed",frk->name);
    *status=WMAP_TIMEOUT;
    *retry=0;
    *er = addError(wmap_timeout,buferr,*er,__LINE__);
    return;
  }
  return;
}

/* ============================================================ *
 * Translates the input parameter 'parin' to CAMB input params  *
 * 'parout'.							*
 * ============================================================ */
void wmap_reparametrize(common_like *like, const double *parin, double *parout, error** err)
{
   int i;
   double h2;
   double Omegam, Omegab, Omegade, Omeganumass, Omegac, OmegaK,
     omegam, omegab, omegade, omeganumass, omegac, omegaK, h100;
   int iOmegade, iOmegaK, iomegade, iomegaK, i_poly_de;
   wmap_state *state;
   cosmo *model;

   state      = (wmap_state*)(like->state);
   model      = copy_parameters_only(state->model, err);
   forwardError(*err, __LINE__,);

   /* === Default parameters === */
   h2                     = dsqr(model->h_100);
   parout[ci_ns]          = model->n_spec;				/* Scalar index n_s			   */
   /* Delta^2_R has to be set later! */
   parout[ci_Delta2R]     = 0.0;					/* Normalisation Delta^2_R		   */
   parout[ci_alphas]      = 0.0;					/* Running scalar index alpha_s            */
   parout[ci_w0]          = model->w0_de;				/* w_0					   */
   parout[ci_tau]         = 0.087;					/* tau					   */
   parout[ci_Neffnu0]     = 3.04;					/* Number of massless neutrinos Neffnu0    */
   parout[ci_Neffnumass]  = model->Neff_nu_mass;			/* Number of massive neutrinos Neff_numass */
   parout[ci_omeganumass] = model->Omega_nu_mass * h2;			/* omega_nu_mass			   */
   parout[ci_nt]          = 0.0;					/* Spectral index tensors n_t		   */
   parout[ci_r]           = 0.0;					/* Tensor to scalar ratio r		   */
   parout[ci_kpivs]       = 0.002;					/* k pivot scalar			   */
   parout[ci_kpivt]       = 0.002;	                                /* k pivot tensor			   */

   /* Default normalization for extra templates (I expect that the first one is ALWAYS SZ) */  
   for(i=0;i<state->n_template;i++) {
      parout[16+i] = 0;
   }
   state->is_tau = 1;

   reset_base_parameters(&Omegam, &Omegab, &Omegade, &Omeganumass, &Omegac, &OmegaK,
			 &omegam, &omegab, &omegade, &omeganumass, &omegac, &omegaK,
			 &h100, &iOmegade, &iOmegaK, &iomegade, &iomegaK);
   i_poly_de = 0;

   for (i=0; i<like->npar; i++) {
      switch (like->par[i]) {
	 case p_Omegam      : Omegam  = parin[i]; break;
	 case p_Omegab      : Omegab  = parin[i]; break;
	 case p_Omegade     : Omegade = parin[i]; iOmegade = 1; break;
	 case p_Omeganumass : Omeganumass = parin[i]; break;
	 case p_Omegac      : Omegac  = parin[i]; break;
	 case p_OmegaK      : OmegaK  = parin[i]; iOmegaK = 1; break;

	 case p_omegam      : omegam  = parin[i]; break;
	 case p_omegab      : omegab  = parin[i]; break;
	 case p_100_omegab  : omegab  = parin[i]/100.0; break;
	 case p_omegade     : omegade = parin[i]; iomegade = 1; break;
	 case p_omeganumass : omeganumass = parin[i]; break;
	 case p_omegac      : omegac  = parin[i]; break;
	 case p_omegaK      : omegaK  = parin[i]; iomegaK = 1; break;

	 case p_w0de        :
	    testErrorRet(model->de_param == poly_DE, ce_de,
			 "For de_param = 'poly_de', 'w0_de' is not a valid parameter", *err, __LINE__,);
	    parout[ci_w0] = parin[i];
	    break;
	 case p_w1de        : *err = addError(mk_w1de,
					      "Only a constant w_de for CMB is  defined with the currently used version of CAMB",
					      *err, __LINE__);
	    return;
	 case p_wpolyde     :
	    testErrorRetVA(model->N_poly_de <= i_poly_de, ce_de,
			   "N_poly_de (%d) too small", *err, __LINE__,, model->N_poly_de);
	    testErrorRet(i_poly_de > 0, mk_w1de, "Only a constant w_de for CMB is  defined with the currently used version of CAMB", 
			 *err, __LINE__,);
	    model->w_poly_de[i_poly_de] = parin[i];
	    i_poly_de ++;
	    break;

	 case p_h100        : h100 = parin[i]; break; //parout[3] = parin[i]*100.0; break; /* parout[3] is 100*h */
	 case p_sigma8      : *err = addError(mk_sigma8,
					      "Parameter sigma8 for CMB not defined. Instead, use Delta2R "
					      "and add sigma8 as deduced parameter\n", *err, __LINE__);
	    return;
	 case p_Delta2R     : parout[ci_Delta2R]    = parin[i]*1.0e-9; break;
	 case p_Neffnumass  : parout[ci_Neffnumass] = parin[i]; break;
	 case p_Neffnu0     : parout[ci_Neffnu0]    = parin[i]; break;
	 case p_ns          : parout[ci_ns]         = parin[i]; break;
	 case p_tau         : parout[ci_tau]        = parin[i]; break;
	 case p_alphas      : parout[ci_alphas]     = parin[i]; break;
	 case p_nt          : parout[ci_nt]         = parin[i]; break;
	 case p_r           : parout[ci_r]          = parin[i]; break;
	 case p_lnr         : parout[ci_r]          = exp(parin[i]); break;
	 case p_A_SZ        : parout[16]            = parin[i]; break;     /* First template parameter is A_SZ */
	 default            : break;
      }
   }

   if (h100<0) h100 = model->h_100;  /* Default value */
   else model->h_100 = h100;         /* Input parameter value */

   set_base_parameters(model, Omegam, Omegab, Omegade, Omeganumass,
		       Omegac, OmegaK, omegam, omegab, omegade, omeganumass,
		       omegac, omegaK, h100, iOmegade, iOmegaK, iomegade, iomegaK, err);
   forwardError(*err, __LINE__,);


   /* Set density parameters */
   h2                     = dsqr(model->h_100);
   parout[ci_omegac]      = (model->Omega_m - model->Omega_b) * h2;	  /* omega_c				   */
   parout[ci_omegab]      = model->Omega_b * h2;			  /* omega_b				   */
   /* Bug fix (25/09/2012, thanks to F. Simpson): The following line had an erroneous h2. *
    * In camb, omk is Omega_K and not omega_K!						  */
   parout[ci_omegaK]      = 1.0 - model->Omega_m - model->Omega_nu_mass
			     - model->Omega_de; 			  /* Omega_K				   */
   parout[ci_H0]          = model->h_100 * 100;			 	  /* H_0			       	   */
   parout[ci_omeganumass] = model->Omega_nu_mass * h2;                    /* omega_nu_mas			   */
}

/* ============================================================ *
 * Returns CMB log-likelihood by calling the WMAP likelihood    *
 * library functions.						*
 * ============================================================ */
double likeli_CMB(common_like *like, const double* params, error** err)
{
  char buf[2048*100];
  char out_Cl_str[1024], out_Pk_str[1024], out_lensedCl_str[1024], trans_str[4096], nu_str[1024],
    par_str[2048], reion_str[2048];
  int l_max_scalar, lcur,i,nread,l;
  int nok,errV,status;
  double res;
  double *clT,*clE,*clB,*clTE;
  struct timeval before, after, finally;
  long inter, ntimeout;
  int n;
  double sig8;
  char ws8, lens;
  wmap_state *state;

  
  state = (wmap_state*)(like->state);

#ifdef COMM_DEBUG
  fprintf(stderr, "likeli_CMB starting\n");
#endif

  gettimeofday(&before,NULL);

  /* Reparametrize */
  wmap_reparametrize(like, params, state->parout, err);
  forwardError(*err,__LINE__,-1);

  l_max_scalar = state->lmax;
  clT          = state->cls;
  clE          = state->cls+(state->lmax+1);
  clB          = state->cls+(state->lmax+1)*2;
  clTE         = state->cls+(state->lmax+1)*3;
  
  /* Prepare CAMB parameter file */
  if (state->want_sigma_8!=0) {
    ws8='T';
  } else {
    ws8='F';
  }
  if (state->do_lensing!=0) {
    lens='T';
   } else {
    lens='F';
   }
  
  // transfert part
  sprintf(trans_str,"\
CMB_outputscale         = 7.4311e12\n\
get_scalar_cls          = T\n\
get_tensor_cls          = T\n\
do_lensing              = %c\n\
do_nonlinear            = 0\n\
l_max_scalar            = %d\n\
k_eta_max_scalar        = %d\n\
l_max_tensor            = %d\n\
k_eta_max_tensor        = %d\n\
transfer_kmax           = 0.8\n\
initial_power_num       = 1\n\
initial_condition       = 1\n\
initial_vector          = -1 0 0 0 0\n\
vector_mode             = 0\n\
COBE_normalize          = F\n\
transfer_high_precision = F\n\
transfer_k_per_logint   = 0\n\
transfer_num_redshifts  = 1\n\
transfer_redshift(1)    = 0\n\
scalar_amp(1)           = %g\n\
scalar_spectral_index(1)= %g\n\
scalar_nrun(1)          = %g\n\
tensor_spectral_index(1)= %g\n\
initial_ratio(1)        = %g\n\
pivot_scalar            = %g\n\
pivot_tensor            = %g\n\
get_transfer            = %c\n",
          lens,l_max_scalar,l_max_scalar*2,
          l_max_scalar,(int)(l_max_scalar*5.0/2.0),
          state->parout[5],state->parout[4],
          state->parout[6],state->parout[12],
          state->parout[13],state->parout[14],
          state->parout[15],ws8);
  
  // neutrino part
  sprintf(nu_str,"\
massive_neutrinos       = %g\n\
massless_neutrinos      = %g\n",state->parout[10],state->parout[9]);

  //parameters part
  sprintf(par_str,"\
use_physical            = T\n\
ombh2                   = %g\n\
omch2                   = %g\n\
omnuh2                  = %g\n\
omk                     = %g\n\
hubble                  = %g\n\
w                       = %g\n\
cs2_lam                 = 1\n",
          state->parout[1],state->parout[0],
          state->parout[11],state->parout[2],
          state->parout[3],state->parout[7]);

  
  //Reionization part  
  sprintf(reion_str,"\
reionization            = T\n\
accurate_reionization   = %c\n\
accurate_polarization   = %c\n", TRUE_FALSE(state->accurate), TRUE_FALSE(state->accurate));

  if (state->is_tau!=0) {
    sprintf(reion_str,"%s\n\
re_use_optical_depth    = T\n\
re_optical_depth           = %g\n", reion_str,state->parout[8]);
  } else {    
    sprintf(reion_str,"%s\n\
re_use_optical_depth    = F\n\
re_redshift           = %g\n\
re_delta_redshift    = 1.5\n", reion_str,state->parout[8]);
  }

  /* RECFAST 1.5 recombination parameters */
  sprintf(reion_str, "%s\
RECFAST_fudge = 1.14\n				\
RECFAST_fudge_He = 0.86\n			\
RECFAST_Heswitch = 6\n				\
RECFAST_Hswitch  = T\n", reion_str);

  sprintf(out_lensedCl_str, "%s", "");
  if (state->want_out_Cl==1) {
     sprintf(out_Cl_str, "scalar_output_file      = scalCls.dat\n");
     if (state->do_lensing!=0) {
	sprintf(out_lensedCl_str, "lensed_output_file      = lensedCls.dat\n");
     }
  } else {
     sprintf(out_Cl_str, "%s", "");
  }

  if (state->want_out_Pk==1) {
     sprintf(out_Pk_str, "transfer_matterpower(1) = Pdelta.dat\n");
  } else {
     sprintf(out_Pk_str, "%s", "");
  }
  
  
  /* Combine all strings */
  sprintf(buf, "number_of_threads = 0\n%s\n%s\n%s\n%s\n%s\n%s\n%s\nEND\n",
	  trans_str, nu_str, par_str, reion_str, out_Cl_str, out_lensedCl_str, out_Pk_str);

  
  /* ============================================================ *
   * CAMB							   *
   * ============================================================ */
  
  nok = 2;
  status = WMAP_NOSTATUS;
  while (nok) {
    nok--;
    
#ifdef COMM_DEBUG
    fprintf(stderr, "likeli_CMB: Calling sendAndCheck for camb\n");
#endif

    errV = sendAndCheck(state->camb,buf,strlen(buf),state->camb->timout);

#ifdef COMM_DEBUG
    fprintf(stderr, "sendAndCheck CAMB, errV=%d (%s) (time-out %ld)\n", errV,
	    errV==0?"success":"failure", (long)state->camb->timout);
#endif

    if (errV==0) {
      int smth;
      smth = 0;
      
      if (alive(state->camb)) {

#ifdef COMM_DEBUG
	 fprintf(stderr, "camb is alive\n");
#endif
        /* Jump over extra outputs of the code */
        while(wmap_readline(state->camb->pout,buf)>=0);

        /* Read the Cl's */
        while (wmap_readline(state->camb->pout,buf)>=0) {
	   n = sscanf(buf,"%d %lg %lg %lg %lg",&(lcur),&(clT[smth+2]),
		      &(clE[smth+2]),&(clB[smth+2]),&(clTE[smth+2]));
	   if (n!=5) {
	      fprintf(stderr, "Bad format read from camb (n=%d, expected n=5 columns)\n", n);
	      smth = 0;
	      break;
	   }
	  //fprintf(stderr, "read: %d\n", buf);
	   smth++; 
        }

#ifdef COMM_DEBUG
	fprintf(stderr, "Read from camb: smth=%d\n", smth);
#endif
        
        /* Read sigma_8 if desired */
        if (state->want_sigma_8==1) {
          while (wmap_readline(state->camb->pout,buf)>=0) {
            n = sscanf(buf," %lg", &sig8);

#ifdef COMM_DEBUG
            fprintf(stderr, "### sig8 = %g (n=%d) ###\n", sig8, n);
#endif

          }
        } else {
	   state->sigma_8 = -1.0;
	}
	state->sigma_8 = sig8;

        pipe_cleanup(state->camb->pout);
      }
      
      if (smth==0) {
        fprintf(stderr, "Camb choking with value (no data read), skipping\n");
        nok    = 0;
        status = WMAP_CHOKE;
        *err = addError(wmap_choke, "Camb choking with value (no data read)", *err, __LINE__);
      } else if (smth==state->lmax-1) {
	 status = WMAP_OK;
	 nok    = 0;
      }
    }
    
    dealWithIssue(state->camb,errV,&nok,&status,err);
    if (isError(*err)) {
      fprintf(stderr,"Camb error, status nok=%d. Make sure, 'scamb' is executable and can be found.\n", nok);
      printError(stderr,*err);
    }
    forwardError(*err,__LINE__,status);
  }
  
  if (status != WMAP_OK) {
#ifdef COMM_DEBUG
    fprintf(stderr, "wmap status = %d not ok, returning prematurely\n", status);
#endif
    return status;
  }
  
  /* Add templates to Cl's. For now only on T */
  if (state->n_template>0) {
    for(l=0;l<state->lmax+1;l++) {
      for(i=0;i<state->n_template;i++) {
        clT[l] += state->parout[16+i] * state->template[l+i*(state->lmax+1)];
      }
    }
  }
  
  gettimeofday(&after,NULL);
  inter = ((after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec))/1000;
  
#ifdef COMM_DEBUG
  fprintf(stderr, "camb used %ld ms\n", (long)inter);
#endif

  if (inter>state->camb->timout) {
#ifdef COMM_DEBUG
    fprintf(stderr, "Time-out after camb\n");  
#endif
    *err = addError(wmap_timeout, "Timeout after camb", *err, __LINE__);
    return WMAP_TIMEOUT;
  }
  ntimeout = state->camb->timout - inter;
  
#ifdef COMM_DEBUG
   fprintf(stderr, "Starting WMAP likelihood\n");
#endif

  
  /* ============================================================ *
   * WMAP{3,5,7} Likelihood					   *
   * ============================================================ */
  
  /* Do something wth the bloody cls! */
  nok    = 4;
  status = WMAP_NOSTATUS;

  while (nok) {
    nok--;

#ifdef COMM_DEBUG
    print_parameter(stderr, like->npar, params);
    fprintf(stderr, "Calling sendAndCheck for wmap likelihood\n");
#endif

    errV = sendAndCheck(state->lkl, (void*)state->cls, sizeof(double)*(state->lmax+1)*4, ntimeout);

#ifdef COMM_DEBUG
    fprintf(stderr, "sendAndCheck WMAP, errV=%d (%s) (time-out %ld)\n", errV,
	    errV==0?"success":"failure", (long)state->camb->timout);
#endif

    if (errV==0) {
      int smth;
      smth = 0;

      if (alive(state->camb)) {
	while(wmap_readline(state->lkl->pout,buf)>=0); /* Jump over extra outputs of the code */
        while (wmap_readline(state->lkl->pout,buf)>=0) {
          nread = sscanf(buf,"%lg",&res);
          if (nread!=1) {
            status = WMAP_CHOKE;
            nok    = 0;
            *err = addErrorVA(wmap_choke, "Bad fomat read from wmap (n=%d, expected n=1)", *err, __LINE__, nread);
            break;
          }
          smth = 1;
        }

        pipe_cleanup(state->lkl->pout);
      }

      if (smth==0) {
        *err   = addError(wmap_choke, "wmap_lkl choking with value (no data read)", *err, __LINE__);
        nok    = 0;
        status = WMAP_CHOKE;
      } else if (smth==1) {
	 status = WMAP_OK;
	 nok    = 0;
      }
    }

    dealWithIssue(state->lkl,errV,&nok,&status,err);
    forwardError(*err,__LINE__,status);
  }
  
  gettimeofday(&finally,NULL);
#ifdef COMM_DEBUG
  fprintf(stderr, "time spent: %d msec (camb %d, wmap_lkl %d)\n",
          abs((before.tv_sec-finally.tv_sec)*1000000+(before.tv_usec-finally.tv_usec))/1000,
          abs((after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec))/1000,
          abs((after.tv_sec-finally.tv_sec)*1000000+(after.tv_usec-finally.tv_usec))/1000);
#endif
  
  if (status != WMAP_OK) { 
    fprintf(stderr, "wmap status = %d not ok, returning late\n", status);
    return status;
  }
  purgeError(err);

#ifdef COMM_DEBUG  
  fprintf(stderr, "wmap status = %d, returning res = %g\n", status, res);
#endif
  
  return res;
}

void deduced_CMB(const common_like *like, double *res)
{
   wmap_state *state;

   state = (wmap_state*)like->state;

   res[0] = state->sigma_8;
}

void print_CMB(FILE *where, void *state, error **err)
{
   wmap_state *wmapstate;
   FILE *rhere;

   wmapstate = (wmap_state*)state;

   if (where==NULL) rhere = stdout;
   else rhere = where;

   fprintf(rhere, "scamp_path   = %s\n", wmapstate->scamb_path);
   fprintf(rhere, "data_path    = %s\n", wmapstate->data_path);
   fprintf(rhere, "Cl_SZ_file   = %s\n", wmapstate->Cl_SZ_file);
   fprintf(rhere, "model_file   = %s\n", wmapstate->model_file);
   fprintf(rhere, "want_sigma_8 = %d\n", wmapstate->want_sigma_8);
   fprintf(rhere, "want_out_Cl  = %d\n", wmapstate->want_out_Cl);
   fprintf(rhere, "want_out_Pk  = %d\n", wmapstate->want_out_Pk);
   fprintf(rhere, "do_lensing   = %d\n", wmapstate->do_lensing);
   fprintf(rhere, "is_tau       = %d\n", wmapstate->is_tau);
   fprintf(rhere, "n_template   = %d\n", wmapstate->n_template);
   fprintf(rhere, "lmax         = %d\n", wmapstate->lmax);
   fprintf(rhere, "accurate     = %d\n", wmapstate->accurate);
}

#define NCWDBUF 2048
void check_scamb(void *state, error **err)
{
   char cwd[NCWDBUF];
   FILE *SCAMB;
   wmap_state *wmapstate;

   wmapstate = (wmap_state*)state;

   /* See wether scamb executable exists */
   SCAMB = fopen(wmapstate->scamb_path, "r");
   testErrorRetVA(SCAMB==NULL, mk_scamb,
		  "'%s' does not exist, see 'scamb_path' key in config file",
		  *err, __LINE__,, wmapstate->scamb_path);
   fclose(SCAMB);

   /* Check the WMAP directory */
   testErrorRet(getcwd(cwd, NCWDBUF)==NULL, mk_cwdbuf,
		"Buffer to store cwd() is not large enough", *err, __LINE__,);
   testErrorRetVA(chdir(wmapstate->data_path)!=0, mk_chdir,
		  "Cannot change into WMAP data directory '%s'\n"
		  "(see 'data_path' key in config file, CMB section)",
		  *err, __LINE__,, wmapstate->data_path);
   chdir(cwd);
}


/* ============================================================ *
 * WMAP5 distance priors (Komatsu et al. 2008)			*
 * ============================================================ */

functions_wrapper_t *init_functions_CMBDistPrior(error **err)
{
   functions_wrapper_t *init_func;

   //fprintf(stderr, "init_functions_CMBDistPrior\n");

   init_func = init_functions_wrapper(read_from_config_CMBDistPrior, init_CMBDistPrior, likeli_CMBDistPrior,
				      special_CMBDistPrior, NULL, print_CMBDistPrior, err);
   forwardError(*err, __LINE__, NULL);

   return init_func;
}

void read_from_config_CMBDistPrior(void **state, FILE *F, error **err)
{
   cmbDP_state *cmbDPstate;
   config_element c = {0, 0.0, ""};
   int j;

   cmbDPstate = malloc_err(sizeof(cmbDP_state), err);
   forwardError(*err, __LINE__,);

   CONFIG_READ_S(cmbDPstate, datname, s, F, c, err);
   CONFIG_READ_S(cmbDPstate, model_file, s, F, c, err);

   CONFIG_READ_S(cmbDPstate, sspecial, s, F, c, err);
   STRING2ENUM(cmbDPstate->special, cmbDPstate->sspecial, special_t, sspecial_t, j, Nspecial_t, err);

   *state = cmbDPstate;
}

void init_CMBDistPrior(common_like *like, error **err)
{
   cmbDP_state *cmbDPstate;
   FILE *DATA, *FD;

   cmbDPstate  = (cmbDP_state*)(like->state);

   /* Read data file with ML shift parameters and inverse covariance */
   DATA = fopen_err(cmbDPstate->datname, "r", err);              forwardError(*err, __LINE__,);
   cmbDPstate->data = mvdens_dwnp(DATA, err);           forwardError(*err, __LINE__,);
   fclose(DATA);
   //sm2_inverse(cmbDPstate->data->std, cmbDPstate->data->ndim, err);
   mvdens_inverse(cmbDPstate->data, err);
   forwardError(*err, __LINE__,);

   /* Default cosmological model */
   if (strcmp(cmbDPstate->model_file, "-")!=0) {
      FD = fopen_err(cmbDPstate->model_file, "r", err);
      forwardError(*err, __LINE__,);
      read_cosmological_parameters(&cmbDPstate->model, FD, err);
      forwardError(*err, __LINE__,);
      fclose(FD);
   } else {
      cmbDPstate->model = set_cosmological_parameters_to_default2(err);
      forwardError(*err, __LINE__,);
   }
}

#define EPS              1.0e-6
double likeli_CMBDistPrior(common_like *like, const double *params, error **err)
{
   cmbDP_state *cmbDPstate;
   cosmo *model;
   double res;
   double Omegam, Omegab, Omegade, Omeganumass, Omegac, OmegaK,
     omegam, omegab, omegade, omeganumass, omegac, omegaK, h100;
   int iOmegade, iOmegaK, iomegade, iomegaK;
   int i, de_prior, i_poly_de;

   cmbDPstate = (cmbDP_state*)(like->state);

   model = copy_parameters_only(cmbDPstate->model, err);
   forwardError(*err, __LINE__, 0.0);
   //dump_param(model, stderr);

   reset_base_parameters(&Omegam, &Omegab, &Omegade, &Omeganumass, &Omegac, &OmegaK,
			 &omegam, &omegab, &omegade, &omeganumass, &omegac, &omegaK,
			 &h100, &iOmegade, &iOmegaK, &iomegade, &iomegaK);
   i_poly_de = 0;

   for (i=0; i<like->npar; i++) {
      testErrorRet(!finite(params[i]), ce_infnan, "params[i] not finite", *err, __LINE__, 0.0);
      switch (like->par[i]) {
	 case p_Omegam      : Omegam  = params[i]; break;
	 case p_Omegab      : Omegab  = params[i]; break;
	 case p_Omegade     : Omegade = params[i]; iOmegade = 1; break;
	 case p_Omeganumass : Omeganumass = params[i]; break;
	 case p_Omegac      : Omegac  = params[i]; break;
	 case p_OmegaK      : OmegaK  = params[i]; iOmegaK = 1; break;

	 case p_omegam      : omegam  = params[i]; break;
	 case p_omegab      : omegab  = params[i]; break;
	 case p_100_omegab  : omegab  = params[i]/100.0; break;
	 case p_omegade     : omegade = params[i]; iomegade = 1; break;
	 case p_omeganumass : omeganumass = params[i]; break;
	 case p_omegac      : omegac  = params[i]; break;
	 case p_omegaK      : omegaK  = params[i]; iomegaK = 1; break;

	 case p_w0de        :
	    testErrorRet(model->de_param == poly_DE, ce_de,
			 "For de_param = 'poly_de', 'w0_de' is not a valid parameter", *err, __LINE__, 0.0);
	    model->w0_de = params[i];
	    break;
	 case p_w1de        :
	    testErrorRet(model->de_param == poly_DE, ce_de,
			 "For de_param = 'poly_de', 'w0_de' is not a valid parameter", *err, __LINE__, 0.0);
	    model->w1_de = params[i];
	    break;
	 case p_wpolyde     :
	    testErrorRetVA(model->N_poly_de <= i_poly_de, ce_de,
			   "N_poly_de (%d) too small", *err, __LINE__, 0.0, model->N_poly_de);
	    model->w_poly_de[i_poly_de] = params[i];
	    i_poly_de ++;
	    break;

	 case p_h100        : h100 = params[i]; break;
	 case p_sigma8      : model->sigma_8 = params[i]; break;
	 case p_Delta2R     : *err = addError(mk_sigma8,
					      "Parameter Delta2R for CMBDistPrior not defined. Instead, use sigma8",
					      *err, __LINE__);
	    return 0.0;
	 case p_ns          : model->n_spec = params[i]; break;
	 case p_Neffnumass  : model->Neff_nu_mass = params[i]; break;
	 default            : break;
      }
   }


   if (h100<0) h100 = model->h_100;  /* Default value */
   else model->h_100 = h100;         /* Input parameter value */

   set_base_parameters(model, Omegam, Omegab, Omegade, Omeganumass,
		       Omegac, OmegaK, omegam, omegab, omegade, omeganumass,
		       omegac, omegaK, h100, iOmegade, iOmegaK, iomegade, iomegaK, err);
   forwardError(*err, __LINE__, -1.0);


   /* Delete precomputed tables in model if parameters have changed */
   updateFrom(cmbDPstate->model, model, err);
   forwardError(*err, __LINE__, 0.0);

   de_prior = 0;
   if (cmbDPstate->special==de_conservative) {
      de_prior = test_range_de_conservative(model, err);
      forwardError(*err, __LINE__, 0.0);
   }

   if (de_prior==0) {
      res = chi2_cmbDP(model, cmbDPstate->data, err);
      forwardError(*err, __LINE__, 0.0);
   } else {
      res = 0.0;
   }

   // fprintf(stderr, "res2 = %g\n", res);
   free_parameters(&model);

   if (de_prior==1) {
      *err = addError(wmap_de_prior, "Dark-energy eos w outside of conservative prior", *err, __LINE__);
      return 0.0;
   }

   return res;
}
#undef EPS

special_t special_CMBDistPrior(void *state)
{
   cmbDP_state *cmbDPstate;
   cmbDPstate = (cmbDP_state*)state;
   return cmbDPstate->special;
}

void print_CMBDistPrior(FILE *where, void *state, error **err)
{
   FILE *rhere;
   cmbDP_state *cmbDPstate;

   rhere = where;
   if (where==NULL) {
      rhere = stdin;
   }

   cmbDPstate = (cmbDP_state*)state;

   fprintf(rhere, "datname = %s\n", cmbDPstate->datname);
}


/* ============================================================ *
 * Dummy functions if WMAP is not defined in Makefile.		*
 * The WMAP-counterparts are in CMB/C_wrappers.c		*
 * ============================================================ */

#ifndef DOWMAP

void* lkl_init(const char *dummy, error **err)
{
   fflush(stderr);
   fprintf(stderr, "Error: CMB not supported. Re-run 'configure.py' with option '-c'.\n");
   exit(2);
}
 
double compute_lkl(void* ptr, double* tt, double* te, double* ee, double* bb, 
		   int *pnlmax)
{
   fflush(stderr);
   fprintf(stderr, "Error: CMB not supported. Re-run 'configure.py' with option '-c'.\n");
   exit(2);
}
 
void f90_compute_lkl_(double* cl_tt, double* cl_te, double *cl_ee, double *cl_bb,
		      int *nlmax, double *like_tot)
{}

#endif /* ifndef DOWMAP */


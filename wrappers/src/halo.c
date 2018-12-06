/* ============================================================ *
 * halo.c	                                    					 *
 * Martin Kilbinger, Jean Coupon                                *
 * 2009 - 2014			                                           *
 * ============================================================ */

#include "halo.h"


functions_wrapper_t *init_functions_GalCorr(error **err)
{
   functions_wrapper_t *init_func;

   init_func = init_functions_wrapper(read_from_config_GalCorr, init_GalCorr, likeli_GalCorr,
				      NULL, NULL, print_GalCorr, err);
   forwardError(*err, __LINE__, NULL);

   return init_func;
}

void read_from_config_GalCorr(void **state, FILE *F, error **err)
{
   halo_state *halostate;
   config_element c = {0, 0.0, ""};
   int j;
   
   halostate = malloc_err(sizeof(halo_state),err);
   forwardError(*err,__LINE__,);

   CONFIG_READ_S(halostate, shalodata, s, F, c, err);
   STRING2ENUM(halostate->halodata, halostate->shalodata, halodata_t, shalodata_t, j, Nhalodata_t, err);

   CONFIG_READ_S(halostate, shalomode, s, F, c, err);
   STRING2ENUM(halostate->halomode, halostate->shalomode, halomode_t, shalomode_t, j, Nhalomode_t, err);

   CONFIG_READ_S(halostate, datname, s, F, c, err);

   if (halostate->halomode==galcorr_cov || halostate->halomode==galcorr_log) {
      CONFIG_READ_S(halostate, covname, s, F, c, err);
      CONFIG_READ(halostate, corr_invcov, d, F, c, err);
   } else {
      strcpy(halostate->covname, "");
   }

   /* Integral constraint */
   CONFIG_READ_S(halostate, sintconst_type, s, F, c, err);
   STRING2ENUM(halostate->intconst_type, halostate->sintconst_type, intconst_t, sintconst_t, j, Nintconst_t, err);
   if (halostate->intconst_type == constant) {
     CONFIG_READ(halostate, delta, d, F, c, err);
     CONFIG_READ(halostate, intconst, d, F, c, err);
   }else if(halostate->intconst_type == random_file) {
     CONFIG_READ_S(halostate, intconst_file, s, F, c, err);
   }
   
   /* Ngal */
   CONFIG_READ_S(halostate, sngal_fit_type, s, F, c, err);
   STRING2ENUM(halostate->ngal_fit_type, halostate->sngal_fit_type, ngal_fit_t, sngal_fit_t, j, Nngal_fit_t, err);
   if (halostate->ngal_fit_type!=ngal_no_fit) {
     CONFIG_READ(halostate, ngal, d, F, c, err);
     CONFIG_READ(halostate, ngalerr, d, F, c, err);
   } else {
     halostate->ngal = halostate->ngalerr = -1.0;
   }
   
   CONFIG_READ_S(halostate, model_file, s, F, c, err);

   CONFIG_READ_S(halostate, sspecial, s, F, c, err);
   STRING2ENUM(halostate->special, halostate->sspecial, special_t, sspecial_t, j, Nspecial_t, err);

   *state = (void*)halostate;
}

void init_GalCorr(common_like *like, error **err)
{
   halo_state *halostate;
   FILE *FD;
   datcov *dc;
   int i, j;

   halostate = (halo_state*)(like->state);
   
   switch (halostate->halodata) {
     
   case w_of_theta :  case wp_rp : case deltaSigma : case smf :
     
     /* Read data */
     halostate->wt = read_wtheta(halostate->datname, halostate->intconst_type, halostate->delta, halostate->intconst, halostate->intconst_file, err);
     
     forwardError(*err, __LINE__, );
     
     if (halostate->halomode==galcorr_cov || halostate->halomode==galcorr_log) {
       
       /* Read covariance */
       dc = init_datcov_for_cov_only(1, halostate->wt->nbins, err);    forwardError(*err, __LINE__,);
       read_cov_tomo(dc, halostate->covname, 0, err);                  forwardError(*err, __LINE__,);
       
       /* DEBUGGING: [Jean], this is useless */
       if (halostate->halomode==galcorr_log) {
	 /* Transform to log(w), Cov[log(w)] */
	 for (i=0; i<dc->n; i++) {
	   testErrorRetVA(halostate->wt->w[i]<0, hm_negative, "Negative correlation w(theta_%d)=%g not allowed for 'galcorr_log' halomode",
			  *err, __LINE__,, i, halostate->wt->w[i]);
	   halostate->wt->w[i] = log(halostate->wt->w[i]);
	   for (j=0; j<dc->n; j++) {
	     /* sigma(logxi, logxj) = sigma(xi, xj)/(xi*xj) */
	     dc->cov[0][i*dc->n+j] = dc->cov[0][i*dc->n+j]/(halostate->wt->w[i]*halostate->wt->w[j]);
	   }
	 }
       }
       
       //mvdens_inverse(dc, err); // TODO: datcov with mvdens
       sm2_inverse(dc->cov[0], dc->n, err);
       forwardError(*err, __LINE__,);
       
       multiply_all(dc->cov[0], dc->n*dc->n, halostate->corr_invcov);
       //write_matrix(dc->cov[0], dc->n, "cov_inv.dat", err);         forwardError(*err, __LINE__, );
       
       halostate->wt->wcov = malloc_err(dc->n*dc->n*sizeof(double), err);
       forwardError(*err, __LINE__, );
       memcpy(halostate->wt->wcov, dc->cov[0], dc->n*dc->n*sizeof(double));
       del_data_cov(&dc);
     } else {
       halostate->wt->wcov = NULL;
     }
	 break;
	 
   default :
     *err = addErrorVA(ce_unknown, "Unknown or invalid halodata type %d", *err, __LINE__, halostate->halodata);
     return;
   }
   
   /* Cosmological model */
   if (strcmp(halostate->model_file, "-")!=0) {
     FD = fopen_err(halostate->model_file, "r", err);
     forwardError(*err, __LINE__, );
     read_cosmological_parameters_hm(&halostate->model, FD, err);
     forwardError(*err, __LINE__, );
     fclose(FD);
   } else {
     halostate->model = set_cosmological_parameters_to_default_hm(err);
     forwardError(*err, __LINE__, );
   }
}

void fill_parameters_hm(cosmo_hm *halomodel, const double *params, const par_t *par, int npar, error **err)
{
  int i, i_poly_de;
  double Omegam, Omegab, Omegade, Omeganumass, Omegac, OmegaK,
    omegam, omegab, omegade, omeganumass, omegac, omegaK, h100;
  int iOmegade, iOmegaK, iomegade, iomegaK;
  
  reset_base_parameters(&Omegam, &Omegab, &Omegade, &Omeganumass, &Omegac, &OmegaK,
			&omegam, &omegab, &omegade, &omeganumass, &omegac, &omegaK,
			 &h100, &iOmegade, &iOmegaK, &iomegade, &iomegaK);
   i_poly_de = 0;

   for (i=0; i<npar; i++) {
     switch (par[i]) {
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
       testErrorRet(halomodel->cosmo->de_param == poly_DE, ce_de,
		    "For de_param = 'poly_de', 'w0_de' is not a valid parameter", *err, __LINE__,);
       halomodel->cosmo->w0_de = params[i];
       break;
     case p_w1de        :
       testErrorRet(halomodel->cosmo->de_param == poly_DE, ce_de,
		    "For de_param = 'poly_de', 'w0_de' is not a valid parameter", *err, __LINE__,);
       halomodel->cosmo->w1_de = params[i];
       break;
     case p_wpolyde     :
       testErrorRetVA(halomodel->cosmo->N_poly_de <= i_poly_de, ce_de,
		      "N_poly_de (%d) too small", *err, __LINE__,, halomodel->cosmo->N_poly_de);
       halomodel->cosmo->w_poly_de[i_poly_de] = params[i];
       i_poly_de ++;
       break;
       
     case p_h100        : h100 = params[i]; break;
     case p_sigma8      : halomodel->cosmo->normalization = params[i];
       halomodel->cosmo->normmode = norm_s8;
       break;
     case p_ns          : halomodel->cosmo->n_spec = params[i]; break;
       
     case p_Neffnumass  : halomodel->cosmo->Neff_nu_mass = params[i]; break;
     case p_Mmin        : halomodel->log10M_min  = log10(params[i]);   break;
     case p_log10Mmin   : halomodel->log10M_min  = params[i];          break;
     case p_M1          : halomodel->log10M1     = log10(params[i]);   break;
     case p_log10M1     : halomodel->log10M1     = params[i];          break;
     case p_M0          : halomodel->log10M0     = log10(params[i]);   break;
     case p_log10M0     : halomodel->log10M0     = params[i];          break;
     case p_sigma_log_M : halomodel->sigma_log_M = params[i];          break;
     case p_log10Mhalo  : halomodel->log10Mhalo  = params[i];          break; /* for g-g-lensing */
     case p_alpha_halo  : halomodel->alpha       = params[i];          break;
     case p_alpha       : *err = addError(hm_alpha, "You probably meant 'alpha_halo' instead of 'alpha'"
					  " in the 'spar' parameter list in the config file", *err, __LINE__);
     case p_eta         : halomodel->eta = params[i];                 break;

     /* Additional parameters for leauthaud11 and coupon15 models */
     case p_Mstar0          : halomodel->log10Mstar0 = log10(params[i]);  break;
     case p_log10Mstar0     : halomodel->log10Mstar0 = params[i];         break;
     case p_beta            : halomodel->beta        = params[i];         break;
     case p_delta           : halomodel->delta       = params[i];         break;
     case p_gamma           : halomodel->gamma       = params[i];         break;
     case p_B_cut           : halomodel->B_cut       = params[i];         break;
     case p_B_sat           : halomodel->B_sat       = params[i];         break;
     case p_beta_cut        : halomodel->beta_cut    = params[i];         break;
     case p_beta_sat        : halomodel->beta_sat    = params[i];         break;
     case p_fcen1           : halomodel->fcen1       = params[i];         break;
     case p_fcen2           : halomodel->fcen2       = params[i];         break;

     default            : break;
     }
   }
   
   if (h100<0) h100 = halomodel->cosmo->h_100;  /* Default value */
   else halomodel->cosmo->h_100 = h100;         /* Input parameter value */
   
   set_base_parameters(halomodel->cosmo, Omegam, Omegab, Omegade, Omeganumass,
		       Omegac, OmegaK, omegam, omegab, omegade, omeganumass,
		       omegac, omegaK, h100, iOmegade, iOmegaK, iomegade, iomegaK, err);
   forwardError(*err, __LINE__,);
}

double likeli_GalCorr(common_like *like, const double *params, error **err)
{
   halo_state *halostate;
   cosmo_hm *halomodel;
   double res;
   
   halostate  = (halo_state*)(like->state);
   
   halomodel  = copy_parameters_hm(halostate->model, err);
   forwardError(*err, __LINE__, -1.0);
   
   fill_parameters_hm(halomodel, params, like->par, like->npar, err);
   forwardError(*err, __LINE__, -1.0);
   
   updateFrom_hm(halostate->model, halomodel, err);
   forwardError(*err, __LINE__, -1.0);
   
   print_parameter(stderr, like->npar, params);
   switch (halostate->halodata) {
     
   case w_of_theta : case wp_rp : case deltaSigma : case smf :
     
     res = chi2_hm(halomodel, halostate->halodata, halostate->halomode, halostate->wt, halostate->ngal_fit_type,
		   halostate->ngal, halostate->ngalerr, halostate->intconst_type, err);
     
     forwardError(*err, __LINE__, -1.0);
     break; 
     
   default :
     *err = addErrorVA(ce_unknown, "Unknown or invalid halodata type %d", *err, __LINE__, halostate->halodata);
     return -1.0;
   }
   
   /* Delete previous model */
   free_parameters_hm(&halostate->model);
   
   /* Copy current model to halo structure */
   halostate->model = copy_parameters_hm(halomodel, err);
   forwardError(*err, __LINE__, 0.0);
   
   /* Delete current model */
   free_parameters_hm(&halomodel);
   
   return res;
}

void print_GalCorr(FILE *where, void *state, error **err)
{
   FILE *rhere;
   halo_state *halostate;

   if (where==NULL) rhere = stdin;
   else rhere = where;

   halostate = (halo_state*)state;

   fprintf(rhere, "(s)halodata = (%s)%d\n", halostate->shalodata, halostate->halodata);
   fprintf(rhere, "(s)halomode = (%s)%d\n", halostate->shalomode, halostate->halomode);
   fprintf(rhere, "datname     = %s\n", halostate->datname);
   fprintf(rhere, "covname     = %s\n", halostate->covname);
   fprintf(rhere, "corr_invcov = %g\n", halostate->corr_invcov);
   fprintf(rhere, "delta       = %g\n", halostate->delta);
   fprintf(rhere, "intconst    = %g\n", halostate->intconst);
   /* DEBUGGING */
   //   fprintf(rhere, "area        = %g\n", halostate->area);
   fprintf(rhere, "(s)ngal_fit_type = (%s)%d\n", halostate->sngal_fit_type, halostate->ngal_fit_type);
   fprintf(rhere, "ngal        = %g\n", halostate->ngal);
   fprintf(rhere, "ngalerr     = %g\n", halostate->ngalerr);
   fprintf(rhere, "model_file  = %s\n", halostate->model_file);
   fprintf(rhere, "model:\n");
   dump_param_hm(halostate->model, rhere, err);
   forwardError(*err, __LINE__,);
}

void halo_write_to_config(FILE *OUT, const void *what, error **err)
{
   halo_state *halostate;

   halostate = (halo_state*)(what);

   fprintf(OUT, "# Galaxy clustering data\n");
   fprintf(OUT, "shalodata     %s\n", halostate->shalodata);
   fprintf(OUT, "shalomode     %s\n", halostate->shalomode);
   fprintf(OUT, "datname       %s\n", halostate->datname);
   if (halostate->halomode==galcorr_cov || halostate->halomode==galcorr_log) {
      fprintf(OUT, "covname       %s\n", halostate->covname);
      fprintf(OUT, "corr_invcov   %g\n", halostate->corr_invcov);
   }
   fprintf(OUT, "delta         %g\n", halostate->delta);
   fprintf(OUT, "intconst      %g\n", halostate->intconst);
   
   /* DEBUGGING */
   //   fprintf(OUT, "area          %g   # in degrees\n", halostate->area);
   fprintf(OUT, "(s)ngal_fit_type %s(%d)\n", sngal_fit_t(halostate->ngal_fit_type), halostate->ngal_fit_type);
   fprintf(OUT, "ngal          %g\n", halostate->ngal);
   fprintf(OUT, "ngalerr       %g\n", halostate->ngalerr);
   fprintf(OUT, "(s)special    (%s)%d\n", sspecial_t(halostate->special), halostate->special);
   fprintf(OUT, "model_file    %s\n", halostate->model_file);
}


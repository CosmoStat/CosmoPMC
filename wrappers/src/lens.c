/* ============================================================ *
 * lens.c							*
 * Karim Benabed, Martin Kilbinger				*
 * 2007, 2008							*
 * ============================================================ */

#include "lens.h"

functions_wrapper_t *init_functions_Lensing(error **err)
{
   functions_wrapper_t *init_func;

   init_func = init_functions_wrapper(read_from_config_Lensing, init_Lensing, likeli_Lensing,
				      special_Lensing, NULL, print_Lensing, err);
   forwardError(*err, __LINE__, NULL);

   return init_func;
}

void read_from_config_Lensing(void **state, FILE *F, error **err)
{
   lens_state *lensstate;
   config_element c = {0, 0.0, ""};
   int j;


   lensstate = malloc_err(sizeof(lens_state), err);   forwardError(*err, __LINE__,);

   CONFIG_READ_S(lensstate, slensdata, s, F, c, err);
   STRING2ENUM(lensstate->lensdata, lensstate->slensdata, lensdata_t, slensdata_t, j, Nlensdata_t, err);
   switch (lensstate->lensdata) {
      case xipm: case xip: case xim: case map2poly: case map2gauss: case gsqr: case decomp_eb: case nofz:
         lensstate->order = second_order;
         break;
      case map3gauss: case map3gauss_diag: case map2gauss_map3gauss_diag: case map2gauss_map3gauss:
      case decomp_eb_map3gauss_diag: case decomp_eb_map3gauss:
         lensstate->order = third_order;
         break;
      default:
         *err = addErrorVA(ce_unknown, "Unknown or invalid lensdata type %d(%s)",
               *err, __LINE__, lensstate->lensdata, slensdata_t(lensstate->lensdata));
         return;
   }

   if (lensstate->lensdata==decomp_eb || lensstate->lensdata==decomp_eb_map3gauss_diag ||
         lensstate->lensdata==decomp_eb_map3gauss) {
      CONFIG_READ_S(lensstate, sdecomp_eb_filter, s, F, c, err);
      STRING2ENUM(lensstate->decomp_eb_filter, lensstate->sdecomp_eb_filter, decomp_eb_filter_t,
            sdecomp_eb_filter_t, j, Ndecomp_eb_filter_t, err);

      if (lensstate->decomp_eb_filter == COSEBIs_log) {
         lensstate->cosebi_info = malloc_err(sizeof(cosebi_info_t), err);
         forwardError(*err, __LINE__,);
         CONFIG_READ(lensstate->cosebi_info, th_min, d, F, c, err);
         CONFIG_READ(lensstate->cosebi_info, th_max, d, F, c, err);
         lensstate->cosebi_info->th_min *= arcmin;
         lensstate->cosebi_info->th_max *= arcmin;
         CONFIG_READ_S(lensstate->cosebi_info, path, s, F, c, err);
         testErrorRetVA(is_directory(lensstate->cosebi_info->path) == 0, io_directory,
                        "Path '%s' not a valid directory", *err, __LINE__,, lensstate->cosebi_info->path);
      } else {
         lensstate->cosebi_info = NULL;
      }

   } else {
      lensstate->decomp_eb_filter = decomp_eb_none;
   }

   CONFIG_READ_S(lensstate, sformat, s, F, c, err);
   STRING2ENUM(lensstate->format, lensstate->sformat, lensformat_t, slensformat_t, j, Nlensformat_t, err);
   if (lensstate->format==angle_wquadr) {
      CONFIG_READ(lensstate, a1, d, F, c, err);
      CONFIG_READ(lensstate, a2, d, F, c, err);
   } else {
      lensstate->a1 = lensstate->a2 = 0.0;
   }

   CONFIG_READ_S(lensstate, datname, s, F, c, err);
   if (lensstate->lensdata == map2gauss_map3gauss_diag || lensstate->lensdata == map2gauss_map3gauss
         || lensstate->lensdata == decomp_eb_map3gauss_diag || lensstate->lensdata == decomp_eb_map3gauss) {
      CONFIG_READ_S(lensstate, datname2, s, F, c, err);
   } else {
      strcpy(lensstate->datname2, "-");
   }

   CONFIG_READ_S(lensstate, scov_scaling, s, F, c, err);
   STRING2ENUM(lensstate->cov_scaling, lensstate->scov_scaling, cov_scaling_t, scov_scaling_t, j,
         Ncov_scaling_t, err);

   CONFIG_READ_S(lensstate, covname, s, F, c, err);
   if (lensstate->cov_scaling == cov_ESH09) {

      testErrorRet(lensstate->lensdata != xipm && lensstate->lensdata != xip && lensstate->lensdata != xim,
            lens_cov_scaling,
            "scov_scaling = cov_EHS09 only valid for xi+/xi- covariance", *err, __LINE__,);

     CONFIG_READ_S(lensstate, covname_M, s, F, c, err);
     CONFIG_READ_S(lensstate, covname_D, s, F, c, err);
   } else {
     strcpy(lensstate->covname_D, "-");
     strcpy(lensstate->covname_M, "-");
   }
   /* Set pointers */
   lensstate->covname_ptr[0] = lensstate->covname;
   lensstate->covname_ptr[1] = lensstate->covname_M;
   lensstate->covname_ptr[2] = lensstate->covname_D;

   CONFIG_READ(lensstate, corr_invcov, d, F, c, err);
   testErrorRet(lensstate->cov_scaling == cov_ESH09 && fabs(lensstate->corr_invcov-1)>EPSILON, lens_invcov,
		"For cov_scaling = cov_ESH09, invcov != 1 does not make sense", *err, __LINE__,);

   CONFIG_READ_S(lensstate, model_file, s, F, c, err);

   CONFIG_READ_S(lensstate, sspecial, s, F, c, err);
   STRING2ENUM(lensstate->special, lensstate->sspecial, special_t, sspecial_t, j, Nspecial_t, err);

   *state = (void*)lensstate;
}

void init_Lensing(common_like *like, error** err)
{
   lens_state *lensstate;
   FILE *FD;


   lensstate = (lens_state*)(like->state);

   lensstate->data  = init_data_cov_tomo(lensstate->datname, lensstate->datname2, lensstate->covname_ptr,
					 lensstate->lensdata, lensstate->decomp_eb_filter,
					 lensstate->format, lensstate->corr_invcov,
					 lensstate->a1, lensstate->a2,
					 lensstate->order, lensstate->cov_scaling, err);
   forwardError(*err, __LINE__,);

   lensstate->sigma_8 = -1.0;

   /* Cosmological model */
   if (lensstate->order == second_order) {
      if (strcmp(lensstate->model_file, "-")!=0) {
         FD = fopen_err(lensstate->model_file, "r", err);
         forwardError(*err, __LINE__,);
         read_cosmological_parameters_lens(&lensstate->model, FD, err);
         forwardError(*err, __LINE__,);
         fclose(FD);
      } else {
         lensstate->model = set_cosmological_parameters_to_default_lens(err);
         forwardError(*err, __LINE__,);
      }
      lensstate->model_3rd = NULL;
   } else {
      if (strcmp(lensstate->model_file, "-")!=0) {
         FD = fopen_err(lensstate->model_file, "r", err);
         forwardError(*err, __LINE__,);
         read_cosmological_parameters_lens_3rd(&lensstate->model_3rd, FD, err);
         forwardError(*err, __LINE__,);
         fclose(FD);
      } else {
         lensstate->model_3rd = set_cosmological_parameters_to_default_lens_3rd(err);
         forwardError(*err, __LINE__,);
      }
      lensstate->model = NULL;
   }

   if (lensstate->cov_scaling == cov_ESH09) {
      lensstate->data->fiducial = lensstate->model;
   } else {
      lensstate->data->fiducial = NULL;
   }

   if (lensstate->decomp_eb_filter == COSEBIs_log) {
      lensstate->cosebi_info->n_max = lensstate->data->Ntheta;
   }

}

void fill_parameters_lens(cosmo_lens *lensmodel, cosmo_3rd *lensmodel_3rd, lens_state *lensstate, const double *params, const par_t *par, int npar, error **err)
{
   double Omegam, Omegab, Omegade, Omeganumass, Omegac, OmegaK,
          omegam, omegab, omegade, omeganumass, omegac, omegaK, h100;
   int iOmegade, iOmegaK, iomegade, iomegaK;
   double a_ymmk, b_ymmk, c_ymmk, c0_ymmk0const, a_jonben,
          b_jonben, c_jonben, alpha_ludo, beta_ludo, z0_ludo, z_rescale;
   int i_poly_de, i;

   reset_base_parameters(&Omegam, &Omegab, &Omegade, &Omeganumass, &Omegac, &OmegaK,
			 &omegam, &omegab, &omegade, &omeganumass, &omegac, &omegaK,
			 &h100, &iOmegade, &iOmegaK, &iomegade, &iomegaK);
   a_ymmk =  b_ymmk = c_ymmk = c0_ymmk0const = a_jonben = b_jonben = c_jonben = alpha_ludo = beta_ludo 
     = z0_ludo = z_rescale = -1.0;
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

         case p_w0de        : testErrorRet(lensmodel->cosmo->de_param == poly_DE, ce_de,
                                    "For de_param = 'poly_de', 'w0_de' is not a valid parameter", *err, __LINE__,);
                              lensmodel->cosmo->w0_de = params[i];
                              break;
         case p_w1de        : testErrorRet(lensmodel->cosmo->de_param == poly_DE, ce_de,
                                    "For de_param = 'poly_de', 'w0_de' is not a valid parameter", *err, __LINE__,);
                              lensmodel->cosmo->w1_de = params[i];
                              break;
         case p_wpolyde     : testErrorRetVA(lensmodel->cosmo->N_poly_de <= i_poly_de, ce_de,
                                    "N_poly_de (%d) too small", *err, __LINE__,, lensmodel->cosmo->N_poly_de);
                              lensmodel->cosmo->w_poly_de[i_poly_de] = params[i];
                              i_poly_de ++;
                              break;

         case p_h100        : h100 = params[i]; break;
         case p_sigma8      : lensmodel->cosmo->normalization = params[i];
                              lensmodel->cosmo->normmode = norm_s8;
                              break;
         case p_Delta2R     : testErrorRet(lensstate->sigma_8<0, lens_sigma8, 
                                    "sigma_8 has not yet been calculated (by CMB likelihood)",
                                    *err, __LINE__,);
                              lensmodel->cosmo->normalization = lensstate->sigma_8;
                              /* CAMB has better been called before to set sigma_8! */
                              lensmodel->cosmo->normmode = norm_s8;
                              break;
         case p_ns          : lensmodel->cosmo->n_spec = params[i]; break;
         case p_alphas      : *err = addErrorVA(mk_undef, "Running spectral index '%s' not (yet) implemented for lensing",
                                    *err, __LINE__, spar_t(p_alphas));
                              return;
         case p_Neffnumass  : lensmodel->cosmo->Neff_nu_mass = params[i]; break;
         case p_a_ymmk      : a_ymmk = params[i]; break;
         case p_b_ymmk      : b_ymmk = params[i]; break;
         case p_c_ymmk      : c_ymmk = params[i]; break;
         case p_c0_ymmk0const : c0_ymmk0const = params[i]; break;
         case p_a_jonben    : a_jonben = params[i]; break;
         case p_b_jonben    : b_jonben = params[i]; break;
         case p_c_jonben    : c_jonben = params[i]; break;
         case p_alpha_ludo  : alpha_ludo = params[i]; break;
         case p_beta_ludo   : beta_ludo = params[i]; break;
         case p_z0_ludo     : z0_ludo = params[i]; break;
         case p_z_rescale   : z_rescale = params[i]; break;
         case p_A_ia        :
                              testErrorRet(lensmodel->ia == ia_none, lensing_ia,
                                    "Intrinsic alignment type 'none', but A_ia parameter requested",
                                    *err, __LINE__,);
                              lensmodel->A_ia = params[i];
                              break;

                              /* 3rd-order-only parameters */
         case p_A_GGI       :
                              testErrorRet(lensmodel_3rd->ia == ia_3rd_none, lensing_ia,
                                    "Intrinsic alignment type 'none', but A_GGI parameter requested",
                                    *err, __LINE__,);
                              lensmodel_3rd->A_GGI = params[i];
                              break;
         case p_A_GII       :
                              testErrorRet(lensmodel_3rd->ia == ia_3rd_none, lensing_ia,
                                    "Intrinsic alignment type 'none', but A_GII parameter requested",
                                    *err, __LINE__,);
                              lensmodel_3rd->A_GII = params[i];
                              break;
         case p_theta_GGI       :
                              testErrorRet(lensmodel_3rd->ia == ia_3rd_none, lensing_ia,
                                    "Intrinsic alignment type 'none', but theta_GGI parameter requested",
                                    *err, __LINE__,);
                              lensmodel_3rd->theta_GGI = params[i];
                              break;
         case p_theta_GII       :
                              testErrorRet(lensmodel_3rd->ia == ia_3rd_none, lensing_ia,
                                    "Intrinsic alignment type 'none', but theta_GII parameter requested",
                                    *err, __LINE__,);
                              lensmodel_3rd->theta_GII = params[i];
                              break;
         case p_b_slc :
                              testErrorRet(lensmodel_3rd->slc == slc_none, lensing_3rd_slc,
                                    "Source-lens clustering type sslc='none', but b_slc parameter requested",
                                    *err, __LINE__,);
                              lensmodel_3rd->b_slc = params[i];
	    break;
	 case p_gamma_slc :
	    testErrorRet(lensmodel_3rd->slc == slc_none, lensing_3rd_slc,
			 "Source-lens clustering type sslc='none', but gamma_slc parameter requested",
			 *err, __LINE__,);
	    lensmodel_3rd->gamma_slc = params[i];
	    break;
	 default            : break;
      }
   }


   if (z_rescale>0) {
      /* For the moment: rescale last z-bin */
      lensmodel->redshift->z_rescale[lensmodel->redshift->Nzbin-1] = z_rescale;
   }


   if (h100<0) h100 = lensmodel->cosmo->h_100;  /* Default value */
   else lensmodel->cosmo->h_100 = h100;         /* Input parameter value */

   set_base_parameters(lensmodel->cosmo, Omegam, Omegab, Omegade, Omeganumass,
		       Omegac, OmegaK, omegam, omegab, omegade, omeganumass,
		       omegac, omegaK, h100, iOmegade, iOmegaK, iomegade, iomegaK, err);
   forwardError(*err, __LINE__,);
}


double likeli_Lensing(common_like *like, const double *params, error **err)
{
   lens_state *lensstate;
   cosmo_lens *lensmodel;
   cosmo_3rd *lensmodel_3rd;
   double res;
   int de_prior;


   lensstate  = (lens_state*)(like->state);

   if (lensstate->order == second_order) {
      lensmodel_3rd = NULL;
      lensmodel     = copy_parameters_lens_only(lensstate->model, err);        forwardError(*err, __LINE__, -1.0);
   } else {
      lensmodel_3rd = copy_parameters_3rd_only(lensstate->model_3rd, err);     forwardError(*err, __LINE__, -1.0);
      lensmodel     = lensmodel_3rd->lens;
   }

   fill_parameters_lens(lensmodel, lensmodel_3rd, lensstate, params, like->par, like->npar, err);
   forwardError(*err, __LINE__, -1.0);
   //dump_param_lens(lensmodel, stderr, 0, err); forwardError(*err, __LINE__, -1.0);

   /* Update cosmo parameter in halomodel sub-struct (only relevant if nonlinear=3 (halomodel) */
   copy_parameters_lenshm_cosmo(lensmodel, err);                    forwardError(*err, __LINE__, -1.0);

   /* Delete precomputed tables in lensmodel if parameters in model are different */
   if (lensstate->order == second_order) {
      updateFrom_lens(lensstate->model, lensmodel, err);               forwardError(*err, __LINE__, -1.0);
   } else {
      updateFrom_3rd(lensstate->model_3rd, lensmodel_3rd, err);        forwardError(*err, __LINE__, -1.0);
   }

   de_prior = 0;
   if (lensstate->special==de_conservative) {
      de_prior = test_range_de_conservative(lensmodel->cosmo, err);
      forwardError(*err, __LINE__, -1.0);
   }

   if (de_prior==0) {

      if (lensstate->order == second_order) {
         res = chi2_lensing(lensmodel, lensstate->data, like->want_model, &(like->model), &like->Nmodel,
               lensstate->cosebi_info, err);
         forwardError(*err, __LINE__, -1.0);
      } else {
         res = chi2_lensing_3rd(lensmodel_3rd, lensstate->data, lensstate->cosebi_info, err);
         forwardError(*err, __LINE__, -1.0);
      }

   } else {

      res = 0.0;

   }

   if (lensstate->order == second_order) {
      free_parameters_lens(&lensmodel);
   } else {
      free_parameters_3rd(&lensmodel_3rd);
   }

   return res;
}

special_t special_Lensing(void *state)
{
   lens_state *lensstate;
   lensstate = (lens_state*)state;
   return lensstate->special;
}

void print_Lensing(FILE *where, void *state, error **err)
{
   FILE *rhere;
   lens_state *lensstate;

   if (where==NULL) rhere = stdin;
   else rhere = where;

   lensstate = (lens_state*)state;

   fprintf(rhere, "(s)lensdata         = (%s)%d\n", lensstate->slensdata, lensstate->lensdata);
   fprintf(rhere, "order               = %d\n", lensstate->order);
   fprintf(rhere, "(s)decomp_eb_filter = (%s)%d\n", sdecomp_eb_filter_t(lensstate->decomp_eb_filter),
	   lensstate->decomp_eb_filter);
   if (lensstate->decomp_eb_filter == COSEBIs_log) {
      fprintf(rhere, "th_min              = %g\n", lensstate->cosebi_info->th_min);
      fprintf(rhere, "th_max              = %g\n", lensstate->cosebi_info->th_max);
      fprintf(rhere, "n_max               = %d\n", lensstate->cosebi_info->n_max);
      fprintf(rhere, "path                = %s\n", lensstate->cosebi_info->path);
   }
   fprintf(rhere, "datname             = %s\n", lensstate->datname);
   fprintf(rhere, "datname2            = %s\n", lensstate->datname2);
   fprintf(rhere, "(s)cov_scaling      = (%s)%d\n", lensstate->scov_scaling, lensstate->cov_scaling);
   fprintf(rhere, "covname             = %s\n", lensstate->covname_ptr[0]);
   fprintf(rhere, "covname_D           = %s\n", lensstate->covname_ptr[1]);
   fprintf(rhere, "covname_M           = %s\n", lensstate->covname_ptr[2]);
   fprintf(rhere, "(s)special          = (%s)%d\n", sspecial_t(lensstate->special), lensstate->special);
   fprintf(rhere, "model_file          = %s\n", lensstate->model_file);
   fprintf(rhere, "model:\n");
   if (lensstate->order == second_order) {
      dump_param_lens(lensstate->model, rhere, 1, err);
   } else {
      dump_param_3rd(lensstate->model_3rd, rhere, err);
   }
   forwardError(*err, __LINE__,);
}

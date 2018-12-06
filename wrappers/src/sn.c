/* ============================================================ *
 * sn.c								*
 * Martin Kilbinger 2007					*
 * ============================================================ */

#include "sn.h"

functions_wrapper_t *init_functions_SNIa(error **err)
{
   functions_wrapper_t *init_func;

   //fprintf(stderr, "init_functions_SNIa\n");

   init_func = init_functions_wrapper(read_from_config_SNIa, init_SNIa, likeli_SNIa, special_SN, NULL, print_SNIa, err);
   forwardError(*err, __LINE__, NULL);

   return init_func;
}

void read_from_config_SNIa(void **state, FILE *F, error **err)
{
   Sn_state *snstate;
   config_element c = {0, 0.0, ""};
   int j;
   char s[1024];

   snstate = malloc_err(sizeof(Sn_state), err);   forwardError(*err, __LINE__,);

   CONFIG_READ_S(snstate, datname, s, F, c, err);
   CONFIG_READ_S(snstate, sdatformat, s, F, c, err);
   STRING2ENUM(snstate->datformat, snstate->sdatformat, sndatformat_t, ssndatformat_t, j, Nsndatformat_t, err);

   CONFIG_READ_S(snstate, schi2mode, s, F, c, err);
   STRING2ENUM(snstate->chi2mode, snstate->schi2mode, chi2mode_t, schi2mode_t, j, Nchi2mode_t, err);

   if (snstate->chi2mode==chi2_Theta2_denom_fixed) {
      CONFIG_READ_ARR(snstate, Theta2_denom, d, j, 2, s, F, c, err);
   }

   if (snstate->chi2mode==chi2_dust) {
      CONFIG_READ_S(snstate, zAV_name, s, F, c, err);
      CONFIG_READ_S(snstate, datname_beta_d, s, F, c, err);
   }

   CONFIG_READ(snstate, add_logdetCov, d, F, c, err);

   CONFIG_READ_S(snstate, model_file, s, F, c, err);

   CONFIG_READ_S(snstate, sspecial, s, F, c, err);
   STRING2ENUM(snstate->special, snstate->sspecial, special_t, sspecial_t, j, Nspecial_t, err);

   *state = (void*)snstate;
}

void init_SNIa(common_like *like, error **err)
{
   Sn_state *snstate;
   SnSample *sample;
   FILE *F, *FD;
   int i, n;
   unsigned int Ncomment;
   double z;


   snstate = (Sn_state*)(like->state);

   switch (snstate->datformat) {
      case SNLS_firstyear : case SN_SALT :
	 sample = SnSample_read(snstate->datname, snstate->datformat, err);
	 forwardError(*err, __LINE__, );
	 break;
      default :
	 *err = addErrorVA(sn_unknown, "Unknown datformat %d\n", *err, __LINE__, snstate->datformat);
	 return;
   }
   snstate->sample = sample;


   if (snstate->chi2mode==chi2_dust) {
      n = numberoflines_comments(snstate->zAV_name, &Ncomment, err);
      forwardError(*err, __LINE__, );
      snstate->zAV = init_interTable(n, 0.0, 0.0, 0.0, 0.0, 0.0, err);
      forwardError(*err, __LINE__, );

      F = fopen_err(snstate->zAV_name, "r", err);          forwardError(*err, __LINE__, );
      for (i=0; i<n; i++) {
	 fscanf(F, "%lg %lg\n", &z, snstate->zAV->table+i);
	 if (i==0) snstate->zAV->a   = z;
	 if (i==1) snstate->zAV->dx  = z-snstate->zAV->a;
	 if (i==n-1) snstate->zAV->b = z;
      }
      fclose(F);

      if (strcmp(snstate->datname_beta_d, "-")!=0) {
	 F = fopen_err(snstate->datname_beta_d, "r", err); forwardError(*err, __LINE__, );
	 snstate->data_beta_d = mvdens_dwnp(F, err);       forwardError(*err, __LINE__, );
	 fclose(F);
	 //sm2_inverse(snstate->data_beta_d->std, snstate->data_beta_d->ndim, err);
	 mvdens_inverse(snstate->data_beta_d, err);
	 forwardError(*err, __LINE__, );
      } else {
	 snstate->data_beta_d = NULL;
      }
   } else {
      snstate->zAV         = NULL;
      snstate->data_beta_d = NULL;
   }


   if (snstate->chi2mode==chi2_dust) {
      for (i=0; i<snstate->sample->Nsample; i++) {

	 /* Assign IGM absorption to each object */
	 snstate->sample->data[i].dust = interpol_wr(snstate->zAV, snstate->sample->data[i].z, err);
	 forwardError(*err, __LINE__,);

      }
   } else {
      for (i=0; i<snstate->sample->Nsample; i++) snstate->sample->data[i].dust = 0.0;
   }



   /* Cosmological model */
   if (strcmp(snstate->model_file, "-")!=0) {
      FD = fopen_err(snstate->model_file, "r", err);
      forwardError(*err, __LINE__,);
      read_cosmological_parameters_SN(&snstate->model, FD, err);
      forwardError(*err, __LINE__,);
      fclose(FD);
   } else {
      snstate->model = set_cosmological_parameters_to_default_SN(err);
      forwardError(*err, __LINE__,);
   }
}

#define EPS 1.0e-6
double likeli_SNIa(common_like *like, const double *params, error **err)
{
   Sn_state *snstate;
   cosmo_SN *sn_model;
   double res;
   int i, ibeta_d, iTheta1, de_prior;
   double Omegam, Omegab, Omegade, Omeganumass, Omegac, OmegaK,
     omegam, omegab, omegade, omeganumass, omegac, omegaK,
     beta_d, h100;
   int iOmegade, iOmegaK, iomegade, iomegaK, i_poly_de;


   snstate = (Sn_state*)(like->state);

   //sn_model = set_cosmological_parameters_to_default_SN(err);
   sn_model = copy_parameters_SN_only(snstate->model, err);
   forwardError(*err, __LINE__, 0);

   reset_base_parameters(&Omegam, &Omegab, &Omegade, &Omeganumass, &Omegac, &OmegaK,
			 &omegam, &omegab, &omegade, &omeganumass, &omegac, &omegaK,
			 &h100, &iOmegade, &iOmegaK, &iomegade, &iomegaK);
   beta_d = -1.0;
   OmegaK = omegaK = 0.0;
   ibeta_d = iTheta1 = 0;

   sn_model->stretch = 1.0;
   sn_model->color = 0.0;
   i_poly_de = 0;

   for (i=0; i<like->npar; i++) {
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
	 case p_omegac      : omegac  = params[i]; break;
	 case p_omegaK      : omegaK  = params[i]; iomegaK = 1; break;

	 case p_w0de        :
	    testErrorRet(sn_model->cosmo->de_param == poly_DE, ce_de,
			 "For de_param = 'poly_de', 'w0_de' is not a valid parameter", *err, __LINE__, 0.0);
	    sn_model->cosmo->w0_de = params[i];
	    break;
	 case p_w1de        :
	    testErrorRet(sn_model->cosmo->de_param == poly_DE, ce_de,
			 "For de_param = 'poly_de', 'w0_de' is not a valid parameter", *err, __LINE__, 0.0);
	    sn_model->cosmo->w1_de = params[i];
	    break;
	 case p_wpolyde     :
	    testErrorRetVA(sn_model->cosmo->N_poly_de <= i_poly_de, ce_de,
			   "N_poly_de (%d) too small", *err, __LINE__, 0.0, sn_model->cosmo->N_poly_de);
	    sn_model->cosmo->w_poly_de[i_poly_de] = params[i];
	    i_poly_de ++;
	    break;

	 case p_h100        : h100 = params[i]; break;

	 case p_M           : sn_model->Theta2[0] = params[i]; break;
	 case p_alpha       : sn_model->Theta2[1] = params[i]; break;
	 case p_beta        : sn_model->Theta2[2] = params[i]; break;
	 case p_logbeta     : sn_model->Theta2[2] = -exp(params[i]); break;
	 case p_beta_z      : sn_model->Theta2[3] = params[i]; break;

	 case p_stretch     : sn_model->stretch = params[i]; break;
	 case p_color       : sn_model->color = params[i]; break;

	 case p_theta10     : sn_model->Theta1[0] = params[i]; iTheta1 = 1; break;
	 case p_theta11     : sn_model->Theta1[1] = params[i]; iTheta1 = 1; break;
	 case p_theta12     : sn_model->Theta1[2] = params[i]; iTheta1 = 1; break;
	 case p_theta13     : sn_model->Theta1[3] = params[i]; iTheta1 = 1; break;
	 case p_theta14     : sn_model->Theta1[4] = params[i]; iTheta1 = 1; break;
	 case p_theta15     : sn_model->Theta1[5] = params[i]; iTheta1 = 1; break;
	 case p_theta16     : sn_model->Theta1[6] = params[i]; iTheta1 = 1; break;
	 case p_theta17     : sn_model->Theta1[7] = params[i]; iTheta1 = 1; break;

	 case p_beta_d      : beta_d = params[i]; ibeta_d = 1; break;

	 default : break;
      }
   }

   /* Note that h is absorbed in M, so it's not contrained by SNIa data */

   if (h100<0) h100 = sn_model->cosmo->h_100;  /* Default value */
   else sn_model->cosmo->h_100 = h100;         /* Input parameter value */


   set_base_parameters(sn_model->cosmo, Omegam, Omegab, Omegade, Omeganumass,
		       Omegac, OmegaK, omegam, omegab, omegade, omeganumass,
		       omegac, omegaK, h100, iOmegade, iOmegaK, iomegade, iomegaK, err);
   forwardError(*err, __LINE__, -1.0);


   sn_model->chi2mode        = snstate->chi2mode;
   sn_model->Theta2_denom[1] = snstate->Theta2_denom[1];
   sn_model->Theta2_denom[2] = snstate->Theta2_denom[2];

   if (sn_model->chi2mode==chi2_dust) {
      if (ibeta_d==0) {
	 /* Dust absoption coefficient beta_d not in parameter list: *
	  * Set to best-fit value (e.g. from Menard et al. 2009)     */
	 if (snstate->data_beta_d!=NULL) {
	    sn_model->beta_d = snstate->data_beta_d->mean[0];
	 } else {
	    sn_model->beta_d = 4.9;
	 }
      } else {
	 sn_model->beta_d = beta_d;
      }
   } else {
      /* No dust modelling. Set to zero (to be on safe side) */
      sn_model->beta_d = 0.0;
   }

   /* Pre-calculate luminosity distances for all SNIa */
   SetDl(sn_model, snstate->sample, err);
   forwardError(*err, __LINE__, -1);

   de_prior = 0;
   if (snstate->special==de_conservative) {
      de_prior = test_range_de_conservative(sn_model->cosmo, err);
      forwardError(*err, __LINE__, 0.0);
   }

   if (de_prior==0) {
      res = chi2_SN(sn_model, snstate->sample, snstate->data_beta_d, iTheta1, snstate->add_logdetCov, err);
      forwardError(*err, __LINE__, -1);
   } else {
      res = 0.0;
   }

   //dump_param_SN(sn_model, stderr);  printf("chi2 = %g\n", res); exit(0);

   free_parameters_SN(&sn_model);

   return res;
}
#undef EPS

special_t special_SN(void *state)
{
   Sn_state *snstate;
   snstate = (Sn_state*)state;
   return snstate->special;
}

void print_SNIa(FILE *where, void *state, error **err)
{
   FILE *rhere;
   Sn_state *snstate;

   snstate = (Sn_state*)state;
    
   rhere=where;
   if (where==NULL) {
      rhere=stdin;
   }
   fprintf(rhere, "datname = %s\n", snstate->datname);
   fprintf(rhere, "(s)datformat = (%s)%d\n", snstate->sdatformat, snstate->datformat);
   fprintf(rhere, "(s)chi2mode = (%s)%d\n", snstate->schi2mode, snstate->chi2mode);
   if (snstate->chi2mode==chi2_Theta2_denom_fixed) {
      fprintf(rhere, "Theta2_denom = %f %f\n", snstate->Theta2_denom[0], snstate->Theta2_denom[1]);
   } else if (snstate->chi2mode==chi2_dust) {
      fprintf(rhere, "zAV_name = %s\n", snstate->zAV_name);
      fprintf(rhere, "datname_beta_d = %s\n", snstate->datname_beta_d);
   }
   fprintf(rhere, "add_logdetCov = %d\n", snstate->add_logdetCov);
}




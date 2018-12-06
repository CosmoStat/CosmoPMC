/* ============================================================ *
 * br.c								*
 * Martin Kilbinger						*
 * 2013								*
 * ============================================================ */


#include "br.h"

functions_wrapper_t *init_functions_bias(error **err)
{
   functions_wrapper_t *init_func;

   init_func = init_functions_wrapper(read_from_config_bias, init_bias, likeli_bias,
				      NULL, NULL, print_bias, err);
   forwardError(*err, __LINE__, NULL);

   return init_func;
}

void read_from_config_bias(void **state, FILE *F, error **err)
{
   bias_state *biasstate;
   config_element c = {0, 0.0, ""};
   int j;

   biasstate = malloc_err(sizeof(bias_state), err);
   forwardError(*err, __LINE__,);

   CONFIG_READ_S(biasstate, sdata_type, s, F, c, err);
   STRING2ENUM(biasstate->data_type, biasstate->sdata_type, bias_data_t, sbias_data_t, j, Nbias_data_t, err);

   CONFIG_READ_S(biasstate, datname_gg, s, F, c, err);
   CONFIG_READ_S(biasstate, datname_gn, s, F, c, err);
   CONFIG_READ_S(biasstate, datname_nn, s, F, c, err);
   CONFIG_READ_S(biasstate, covname, s, F, c, err);

   CONFIG_READ_S(biasstate, model_file, s, F, c, err);

   CONFIG_READ_S(biasstate, sspecial, s, F, c, err);
   STRING2ENUM(biasstate->special, biasstate->sspecial, special_t, sspecial_t, j, Nspecial_t, err);


   *state = (void*)biasstate;
}


void init_bias(common_like *like, error **err)
{
   bias_state *biasstate;
   datcov *dc[3];
   int Ntheta_tot, i, k, Nzcorr, Nzbin;
   FILE *FD;

   biasstate = (bias_state*)(like->state);

   biasstate->data = malloc_err(sizeof(bias_data), err);
   biasstate->data->data_type = biasstate->data_type;

   /* Read the three data files */
   for (i=0; i<3; i++) {
      dc[i] = malloc_err(sizeof(datcov), err);
      forwardError(*err, __LINE__,);
   }

   for (i=0; i<3; i++) {
      dc[i]->format = angle_center;
   }

   read_data_tomo(dc[0], biasstate->datname_gg, 0, second_order, err);
   forwardError(*err, __LINE__,);
   biasstate->data->N_gg = dc[0]->Ntheta;
   Nzbin = dc[0]->Nzbin;

   /* Next two call with previous Nzbin, for consistency check */
   read_data_tomo(dc[1], biasstate->datname_gn, Nzbin, second_order, err);
   forwardError(*err, __LINE__,);
   biasstate->data->N_gn = dc[1]->Ntheta;

   read_data_tomo(dc[2], biasstate->datname_nn, Nzbin, second_order, err);
   forwardError(*err, __LINE__,);
   biasstate->data->N_nn = dc[2]->Ntheta;


   /* Copy all three vectors of angular scales */
   biasstate->data->theta_gg = malloc_err(sizeof(double) * biasstate->data->N_gg, err);
   biasstate->data->theta_gn = malloc_err(sizeof(double) * biasstate->data->N_gn, err);
   biasstate->data->theta_nn = malloc_err(sizeof(double) * biasstate->data->N_nn, err);
   for (k=0; k<biasstate->data->N_gg; k++) {
      biasstate->data->theta_gg[k] = dc[0]->theta[k];
   }
   for (k=0; k<biasstate->data->N_gn; k++) {
      biasstate->data->theta_gn[k] = dc[1]->theta[k];
   }
   for (k=0; k<biasstate->data->N_nn; k++) {
      biasstate->data->theta_nn[k] = dc[2]->theta[k];
   }

   /* Copy data vectors into mvdens structure */
   Ntheta_tot             = biasstate->data->N_gg + biasstate->data->N_gn + biasstate->data->N_nn;
   biasstate->data->Nzbin = Nzbin;
   Nzcorr                 = Nzbin * (Nzbin + 1) / 2;
   biasstate->data->mvd = mvdens_alloc(Ntheta_tot * Nzcorr, err);
   forwardError(*err, __LINE__,);

   for (i=k=0; i<biasstate->data->N_gg*Nzcorr; i++,k++) {
      biasstate->data->mvd->mean[k] = dc[0]->data[i];
      //printf("MKDEBUG gg %d %g\n", k, biasstate->data->mvd->mean[k]);
   }
   for (i=0; i<biasstate->data->N_gn*Nzcorr; i++,k++) {
      biasstate->data->mvd->mean[k] = dc[1]->data[i];
      //printf("MKDEBUG gn %d %g\n", k, biasstate->data->mvd->mean[k]);
   }
   for (i=0; i<biasstate->data->N_nn*Nzcorr; i++,k++) {
      biasstate->data->mvd->mean[k] = dc[2]->data[i];
      //printf("MKDEBUG nn %d %g\n", k, biasstate->data->mvd->mean[k]);
   }


   /* Read covariance */

   for (i=0; i<3; i++) del_data_cov(&dc[i]);

   dc[0] = init_datcov_for_cov_only(Nzbin, Ntheta_tot, err);
   forwardError(*err, __LINE__,);
   dc[0]->type   = map2poly; /* Dummy value */


   /* Includes check for consistency with Nzbin, Ntheta_tot from data */
   read_cov_tomo(dc[0], biasstate->covname, 0, err);
   forwardError(*err, __LINE__,);

   // MKDEBUG: inverse, for chi2 by hand
   //fprintf(stderr, "Inverting covariance\n");
   //sm2_inverse(dc[0]->cov[0], dc[0]->n, err);
   //forwardError(*err, __LINE__,);


   for (i=0; i < dc[0]->n * dc[0]->n; i++) {
      biasstate->data->mvd->std[i] = dc[0]->cov[0][i];
   }

   /* Fiducial cosmological model */
   if (strcmp(biasstate->model_file, "-")!=0) {
      FD = fopen_err(biasstate->model_file, "r", err);
      forwardError(*err, __LINE__, );
      read_cosmological_parameters_bias(&biasstate->model, FD, err);
      forwardError(*err, __LINE__, );
      fclose(FD);
   } else {
      biasstate->model = set_cosmological_parameters_to_default_bias(err);
      forwardError(*err, __LINE__, );
   }

   del_data_cov(&dc[0]);
}

void fill_parameters_bias(cosmo_bias *model, const double *params, const par_t *par, int npar, error **err)
{
   int i, i_poly_de, i_kb, i_kr;
   double Omegam, Omegab, Omegade, Omeganumass, Omegac, OmegaK,
     omegam, omegab, omegade, omeganumass, omegac, omegaK, h100;
   int iOmegade, iOmegaK, iomegade, iomegaK;

   reset_base_parameters(&Omegam, &Omegab, &Omegade, &Omeganumass, &Omegac, &OmegaK,
			 &omegam, &omegab, &omegade, &omeganumass, &omegac, &omegaK,
			 &h100, &iOmegade, &iOmegaK, &iomegade, &iomegaK);
   i_poly_de = i_kb = i_kr = 0;

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
	    testErrorRet(model->lens->cosmo->de_param == poly_DE, ce_de,
			 "For de_param = 'poly_de', 'w0_de' is not a valid parameter", *err, __LINE__,);
	    model->lens->cosmo->w0_de = params[i];
	    break;
	 case p_w1de        :
	    testErrorRet(model->lens->cosmo->de_param == poly_DE, ce_de,
			 "For de_param = 'poly_de', 'w0_de' is not a valid parameter", *err, __LINE__,);
	    model->lens->cosmo->w1_de = params[i];
	    break;
	 case p_wpolyde     :
	    testErrorRetVA(model->lens->cosmo->N_poly_de <= i_poly_de, ce_de,
			   "N_poly_de (%d) too small", *err, __LINE__,, model->lens->cosmo->N_poly_de);
	    model->lens->cosmo->w_poly_de[i_poly_de] = params[i];
	    i_poly_de ++;
	    break;

	 case p_h100        : h100 = params[i]; break;
	 case p_sigma8      : model->lens->cosmo->normalization = params[i];
	                      model->lens->cosmo->normmode = norm_s8; break;
	 case p_ns          : model->lens->cosmo->n_spec = params[i]; break;
	 case p_Neffnumass  : model->lens->cosmo->Neff_nu_mass = params[i]; break;

     case p_A_ia        :
        testErrorRet(model->lens->ia == ia_none, lensing_ia,
             "Intrinsic alignment type 'none', but A_ia parameter requested",
             *err, __LINE__,);
        model->lens->A_ia = params[i];
        break;

	 case p_kb :
	    testErrorRetVA(model->Nkb <= i_kb, tls_cosmo_par, // bias_Ncoeff,
			   "Number of bias parameters Nkb=%d too small for %d bias parameters", *err, __LINE__,, model->Nkb, i_kb+1);
	    model->kb[i_kb] = params[i];
	    i_kb ++;
	    break;

	 case p_kr :
	    testErrorRetVA(model->Nkr <= i_kr, tls_cosmo_par,
			   "Number of correlation coefficient parameters Nkr=%d too small for %d correlation parameters",
			   *err, __LINE__,, model->Nkr, i_kr+1);
	    model->kr[i_kr] = params[i];
	    i_kr ++;
	    break;

	 case p_r :
	    *err = addError(hm_alpha, "You probably meant 'kr' instead of 'r'"
			    " in the 'spar' parameter list in the config file",
			    *err, __LINE__);
	    return;

	 default            : break;
      }
   }

   if (h100<0) h100 = model->lens->cosmo->h_100;  /* Default value */
   else model->lens->cosmo->h_100 = h100;         /* Input parameter value */

   set_base_parameters(model->lens->cosmo, Omegam, Omegab, Omegade, Omeganumass,
		       Omegac, OmegaK, omegam, omegab, omegade, omeganumass,
		       omegac, omegaK, h100, iOmegade, iOmegaK, iomegade, iomegaK, err);
   forwardError(*err, __LINE__,);
}

double likeli_bias(common_like *like, const double *params, error **err)
{
   bias_state *biasstate;
   cosmo_bias *model;
   double res;

   biasstate = (bias_state*)(like->state);
   model     = copy_parameters_bias(biasstate->model, err);
   forwardError(*err, __LINE__, -1.0);

   fill_parameters_bias(model, params, like->par, like->npar, err);
   forwardError(*err, __LINE__, -1.0);

   updateFrom_bias(biasstate->model, model, err);
   forwardError(*err, __LINE__, -1.0);

   //print_parameter(stderr, like->npar, params); // MKDEBUG

   res = chi2_bias(model, biasstate->data, err);
   forwardError(*err, __LINE__, -1.0);

   /* Delete previous model */
   free_parameters_bias(&biasstate->model);

   /* Copy current model to bias structure */
   biasstate->model = copy_parameters_bias(model, err);
   forwardError(*err, __LINE__, 0.0);

   /* Delete current model */
   free_parameters_bias(&model);

   return res;
}

void print_bias(FILE *where, void *state, error **err)
{
   FILE *rhere;
   bias_state *biasstate;

   if (where==NULL) rhere = stdin;
   else rhere = where;

   biasstate = (bias_state*)state;

   fprintf(rhere, "(s)dat_atype      = (%s)%d\n", sbias_data_t(biasstate->data_type), biasstate->data_type);
   fprintf(rhere, "datname_gg        = %s\n", biasstate->datname_gg);
   fprintf(rhere, "datname_gn        = %s\n", biasstate->datname_gn);
   fprintf(rhere, "datname_nn        = %s\n", biasstate->datname_nn);
   fprintf(rhere, "covname           = %s\n", biasstate->covname);
   fprintf(rhere, "(s)special        = (%s)%d\n", sspecial_t(biasstate->special), biasstate->special);
   fprintf(rhere, "model_file        = %s\n", biasstate->model_file);
   fprintf(rhere, "model:\n");
   dump_param_bias(biasstate->model, rhere, 1, err);
   forwardError(*err, __LINE__,);
}


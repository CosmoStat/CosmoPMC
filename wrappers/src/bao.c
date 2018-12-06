/* ============================================================ *
 * bao.c							*
 * Martin Kilbinger 2008					*
 * ============================================================ */

#include "bao.h"

functions_wrapper_t *init_functions_BAO(error **err)
{
   functions_wrapper_t *init_func;

   init_func = init_functions_wrapper(read_from_config_BAO, init_BAO, likeli_BAO,
				      special_BAO, NULL, print_BAO, err);
   forwardError(*err, __LINE__, NULL);

   return init_func;
}

void read_from_config_BAO(void **state, FILE *F, error **err)
{
   bao_state *baostate;
   config_element c = {0, 0.0, ""};
   int j;

   baostate = malloc_err(sizeof(bao_state), err);
   forwardError(*err, __LINE__,);

   CONFIG_READ_S(baostate, smethod, s, F, c, err);
   STRING2ENUM(baostate->method, baostate->smethod, method_t, smethod_t, j, Nmethod_t, err);
   CONFIG_READ_S(baostate, datname, s, F, c, err);

   CONFIG_READ_S(baostate, model_file, s, F, c, err);

   CONFIG_READ_S(baostate, sspecial, s, F, c, err);
   STRING2ENUM(baostate->special, baostate->sspecial, special_t, sspecial_t, j, Nspecial_t, err);

   *state = baostate;
}

void init_BAO(common_like *like, error **err)
{
   bao_state *baostate;
   int i, nz;
   FILE *DATA, *FD;
   char str[128];
   config_element c = {0, 0.0, ""};

   baostate = (bao_state*)(like->state);


   /* Read BAO data and redshifts */
   DATA = fopen_err(baostate->datname, "r", err);  forwardError(*err, __LINE__, );
   baostate->data = mvdens_dwnp(DATA, err);        forwardError(*err, __LINE__, );
   nz = baostate->data->ndim;

   if (baostate->method==distance_D_V_ratio) nz = 2*nz;   /* Twice #z as ratios */

   baostate->z = malloc_err(nz*sizeof(double), err); forwardError(*err, __LINE__, );
   CONFIG_READ_ARR(baostate, z, d, i, nz, str, DATA, c, err);
   fclose(DATA);

   mvdens_inverse(baostate->data, err);
   //sm2_inverse(baostate->data->std, baostate->data->ndim, err);
   forwardError(*err, __LINE__, );

   /* Default cosmological model */
   if (strcmp(baostate->model_file, "-")!=0) {
      FD = fopen_err(baostate->model_file, "r", err);
      forwardError(*err, __LINE__, );
      read_cosmological_parameters(&baostate->model, FD, err);
      forwardError(*err, __LINE__, );
      fclose(FD);
   } else {
      baostate->model = set_cosmological_parameters_to_default2(err);
      forwardError(*err, __LINE__, );
   }
} 

#define EPS 1.0e-6
double likeli_BAO(common_like *like, const double *params, error **err)
{
   bao_state *baostate;
   cosmo *model;
   double Omegam, Omegab, Omegade, Omeganumass, Omegac, OmegaK,
     omegam, omegab, omegade, omeganumass, omegac, omegaK, h100;
   double res;
   int iOmegade, iOmegaK, iomegade, iomegaK;
   int i, de_prior, i_poly_de;

   baostate = (bao_state*)(like->state);

   model = copy_parameters_only(baostate->model, err);
   forwardError(*err, __LINE__, -1.0);

   reset_base_parameters(&Omegam, &Omegab, &Omegade, &Omeganumass, &Omegac, &OmegaK,
			 &omegam, &omegab, &omegade, &omeganumass, &omegac, &omegaK,
			 &h100, &iOmegade, &iOmegaK, &iomegade, &iomegaK);
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
   updateFrom(baostate->model, model, err);
   forwardError(*err, __LINE__, 0.0);

   de_prior = 0;
   if (baostate->special==de_conservative) {
      de_prior = test_range_de_conservative(model, err);
      forwardError(*err, __LINE__, 0.0);
   }

   if (de_prior==0) {
      switch (baostate->method) {
	 case distance_A :
	    res = chi2_bao_A(model, baostate->data, baostate->z, err);
	    forwardError(*err, __LINE__, 0.0);
	    break;
	 case distance_d_z :
	    res = chi2_bao_d_z(model, baostate->data, baostate->z, err);
	    forwardError(*err, __LINE__, 0.0);
	    break;
	 case distance_D_V_ratio :
	    res = chi2_bao_D_V_ratio(model, baostate->data, baostate->z, err);
	    break;
	 default :
	    *err = addError(bao_unknown, "Unknown method", *err, __LINE__);
	    return 0.0;
      }
   } else {
      res = 0.0;
   }

   free_parameters(&model);

   return res;
}
#undef EPS

special_t special_BAO(void *state)
{
   bao_state *baostate;
   baostate = (bao_state*)state;
   return baostate->special;
}

void print_BAO(FILE *where, void *state, error **err)
{
   FILE *rhere;
   bao_state *baostate;

   rhere = where;
   if (where==NULL) {
      rhere = stdin;
   }

   baostate = (bao_state*)state;

   fprintf(rhere, "(s)method  = (%s)%d\n", baostate->smethod, baostate->method);
   fprintf(rhere, "datname    = %s\n", baostate->datname);
   fprintf(rhere, "model_file = %s\n", baostate->model_file);
   dump_param(baostate->model, rhere);
}

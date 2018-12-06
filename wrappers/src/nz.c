/* ============================================================ *
 * nz.h								*
 * Martin Kilbinger 2007					*
 * ============================================================ */

#include "nz.h"

functions_wrapper_t *init_functions_Nz(error **err)
{
   functions_wrapper_t *init_func;

   fprintf(stderr, "init_functions_Nz\n");

   init_func = init_functions_wrapper(read_from_config_Nz, init_Nz, likeli_Nz, NULL, NULL, print_Nz, err);
   forwardError(*err, __LINE__, NULL);

   return init_func;
}

void read_from_config_Nz(void **state, FILE *F, error **err)
{
   Nz_state *nzstate;

   nzstate = malloc_err(sizeof(Nz_state), err);
   forwardError(*err, __LINE__,);

   read_redshift_info(&nzstate->redshift, F, err);
   forwardError(*err, __LINE__,);

   *state = nzstate;
}

void init_Nz(common_like *like, error **err)
{
   Nz_state *nzstate;

   nzstate = (Nz_state*)(like->state);

   switch (nzstate->nzmode) {

      case nz_fit_hist :
	 /* Fourth argument = 1: nofz, only variance */
	 // MKDEBUG: TODO, change to init_data_cov_tomo
	 //nzstate->dc = init_data_cov(nzstate->datname, "", "", nofz, 1, err); 
	 nzstate->dc = NULL;
	 forwardError(*err, __LINE__,);
	 break;
      case nz_read_from_files :
	 /* Initialisation of redshift structure was done in Nz_read_from_config */
	 break;
      default :
	 *err = addErrorVA(nz_undef, "Unknown nzmode %d", *err, __LINE__, nzstate->nzmode);
	 return;
   }
}


double likeli_Nz(common_like *like, const double *params, error **err)
{
   redshift_t *redshift;
   Nz_state *nzstate;
   double res, zmin, zmax;
   double a_ymmk, b_ymmk, c_ymmk, a_jonben, b_jonben, c_jonben, alpha_ludo, beta_ludo, z0_ludo;
   int i;

   nzstate = (Nz_state*)(like->state);

   a_ymmk = b_ymmk = c_ymmk = a_jonben = b_jonben = c_jonben = alpha_ludo = beta_ludo = z0_ludo = -1.0;

   for (i=0; i<like->npar; i++) {
      switch (like->par[i]) {
         case p_a_ymmk     : a_ymmk = params[i]; break;
         case p_b_ymmk     : b_ymmk = params[i]; break;
         case p_c_ymmk     : c_ymmk = params[i]; break;
         case p_a_jonben   : a_jonben = params[i]; break;
         case p_b_jonben   : b_jonben = params[i]; break;
         case p_c_jonben   : c_jonben = params[i]; break;
         case p_alpha_ludo : alpha_ludo = params[i]; break;
         case p_beta_ludo  : beta_ludo = params[i]; break;
         case p_z0_ludo    : z0_ludo = params[i]; break;
         default  : break;
      }
   }

   /* TODO: This still only works for one z-bin!! */
   zmin = get_zmin(nzstate->redshift, 0);
   zmax = get_zmax(nzstate->redshift, 0);

   redshift = NULL;

   /* TODO: works only for photo_no */

   /* Redshift distribution */
   if (a_ymmk>0 && b_ymmk>0 && c_ymmk>0) {

      free_and_reinit_redshift(&redshift, 1, 5, err);
      forwardError(*err, __LINE__, 0.0);
      fill_redshift_slice(redshift, 0, ymmk, photz_no, err, zmin, zmax, a_ymmk, b_ymmk, c_ymmk);
      forwardError(*err, __LINE__, 0.0);

   } else if (a_jonben>0 && b_jonben>0 && c_jonben>0) {

      free_and_reinit_redshift(&redshift, 1, 5, err);
      forwardError(*err, __LINE__, 0.0);
      fill_redshift_slice(redshift, 0, jonben, photz_no, err, zmin, zmax, a_jonben, b_jonben, c_jonben);
      forwardError(*err, __LINE__, 0.0);

   } else if (alpha_ludo>0 && beta_ludo>0 && z0_ludo>0) {

      free_and_reinit_redshift(&redshift, 1, 5, err);
      forwardError(*err, __LINE__, 0.0);
      fill_redshift_slice(redshift, 0, ludo, photz_no, err, zmin, zmax, alpha_ludo, beta_ludo, z0_ludo);
      forwardError(*err, __LINE__, 0.0);

   } else {

      *err = addError(lensing_inconsistent, "Unknown n(z) parametrisation",  *err, __LINE__);
      return -1;

   }

   //res = chi2_redshift(redshift, nzstate->dc, err);        forwardError(*err, __LINE__, -1);
   res = 0.0;
   free_redshift(&redshift);

   return res;
}

void print_Nz(FILE *where, void *state, error **err)
{
   FILE *rhere;
   Nz_state *nzstate;

   nzstate = (Nz_state*)state;

   if (where==NULL) rhere = stdout;
   else rhere = where;

   fprintf(rhere, "datname = %s\n", nzstate->datname);
}


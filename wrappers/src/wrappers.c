/* ============================================================ *
 * wrappers.c							*
 * Martin Kilbinger 2010					*		
 * ============================================================ */

#include "wrappers.h"


init_functions_t *init_func_t(data_t data, error **err)
{
   init_functions_t *init_func;

   switch (data) {

      case Mvdens       : init_func = init_functions_Mvdens; break;
      case MixMvdens    : init_func = init_functions_MixMvdens; break;
      case Lensing      : init_func = init_functions_Lensing; break;
      case SNIa         : init_func = init_functions_SNIa; break;
      case Nz           : init_func = init_functions_Nz; break;
      case CMB          : init_func = init_functions_CMB; break;
      case CMBDistPrior : init_func = init_functions_CMBDistPrior; break;
      case BAO          : init_func = init_functions_BAO; break;
	   case Topo         : init_func = init_functions_topo; break;
//    case newmodule    : init_func = init_functions_newmodule; break;
      default           : *err = addErrorVA(wr_undef, "Undefined data type %d(%s)", *err, __LINE__, data, sdata_t(data)); return NULL;

   }

   return init_func;
}


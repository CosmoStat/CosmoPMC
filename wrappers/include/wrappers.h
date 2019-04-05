/* ============================================================ *
 * Martin Kilbinger 2009.					*
 * The mother (file) for all wrappers.				*
 * ============================================================ */

#ifndef __WRAPPERS_H
#define __WRAPPERS_H

#include "pmctools/errorlist.h"
#include "types.h"
#include "init_wrappers.h"
#include "all_wrappers.h"


#define MC_LOGBIG 1.0e10

#define wr_base  -1600
#define wr_undef -1 + wr_base

init_functions_t *init_func_t(data_t data, error **err);

#endif

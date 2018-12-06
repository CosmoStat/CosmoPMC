#ifndef __INIT_WRAPPERS_H
#define __INIT_WRAPPERS_H

#include "errorlist.h"
#include "par.h"

typedef enum {none, unity, de_conservative} special_t;
#define sspecial_t(i) ( \
  i==none ? "none" : \
  i==unity ? "unity" : \
  i==de_conservative ? "de_conservative" :\
  "")
#define Nspecial_t 3

typedef struct {
  par_t *par;
  int npar;
  void *state;
  int want_model, Nmodel;
  double *model;
  struct functions_wrapper_dummy_t *func_wrapper;
} common_like;

typedef void (func_read_from_config_t)(void **, FILE *, error **);
typedef void (func_init_t)(common_like *, error **);
typedef double (func_likeli_t)(common_like *, const double *, error **);
typedef special_t (func_special_t)(void *);
typedef void (func_deduced_t)(const common_like *, double *);
typedef void (func_print_t)(FILE *, void *, error **err);

typedef struct functions_wrapper_dummy_t {
  
  func_read_from_config_t *func_read;
  func_init_t *func_init;
  func_likeli_t *func_likeli;
  func_special_t *func_special;
  func_deduced_t *func_deduced;
  func_print_t *func_print;

} functions_wrapper_t;

typedef functions_wrapper_t *(init_functions_t)(error **);

functions_wrapper_t *init_functions_wrapper(func_read_from_config_t *func_read, func_init_t *func_init,
					    func_likeli_t *func_likeli, func_special_t *func_special,
					    func_deduced_t *func_deduced, func_print_t *func_print, error **err);


#endif

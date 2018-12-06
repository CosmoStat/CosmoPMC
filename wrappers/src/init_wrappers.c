#include "init_wrappers.h"

functions_wrapper_t *init_functions_wrapper(func_read_from_config_t *func_read, func_init_t *func_init,
					    func_likeli_t *func_likeli, func_special_t *func_special, 
					    func_deduced_t *func_deduced, func_print_t *func_print, error **err)
{
   functions_wrapper_t *func_wrapper;

   func_wrapper = malloc_err(sizeof(functions_wrapper_t), err);
   forwardError(*err, __LINE__, NULL);

   func_wrapper->func_read    = func_read;
   func_wrapper->func_init    = func_init;
   func_wrapper->func_likeli  = func_likeli;
   func_wrapper->func_special = func_special;
   func_wrapper->func_deduced = func_deduced;
   func_wrapper->func_print   = func_print;

   return func_wrapper;
}

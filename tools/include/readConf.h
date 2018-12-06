/*
 *  readConf.h
 *  readConf
 *
 *  Created by Karim Benabed on 15/03/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#include "errorlist.h"
#include "io.h"
#include <stdio.h>
#include <string.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include <time.h>


#ifndef __READCONF__
#define __READCONF__
#define RC_BFSTEP 64*1024
#define RC_FNMLEN 4096*4
#define RC_MAX_DEPTH 32
#define RC_MEM_STEP 1024
#define RC_ROOTNAME "rc"
#define RC_DEFAULT_FUNCNAME "rc_dofile"

 
typedef struct {
  void* alias;
  char alias_prefix[RC_FNMLEN+1];
  char fnm[RC_FNMLEN+1];
  void **mem;
  size_t mem_sz;
  size_t mem_cur;
  lua_State *L;
  char last_key[RC_FNMLEN*RC_MAX_DEPTH+RC_MAX_DEPTH+1];
  FILE* display;
  char root[RC_FNMLEN];
  char **read_pars;
} confFile;

confFile* rc_open(char* fnm, error **err);
confFile* rc_init(char* fnm, error **err);
confFile* rc_open_with(char* fnm, char*, int,error **err);

int rc_has_key(confFile* rc, char* key, error **err);

double rc_get_real(confFile *rc, char *key, error **err);
long rc_get_integer(confFile *rc, char *key, error **err);
char* rc_get_string(confFile *rc, char *key, error **err);
double rc_safeget_real(confFile *rc, char *key, double safeguard, error **err);
long rc_safeget_integer(confFile *rc, char *key,long safeguard, error **err);
char* rc_safeget_string(confFile *rc, char *key,char* safeguard, error **err);

int rc_is_array(confFile *rc, char* key, error **err);
size_t rc_array_size(confFile *rc, char* key, error **err);
size_t rc_get_real_array(confFile *rc, char* key, double** parray, error **err);
size_t rc_get_integer_array(confFile *rc, char* key, long** parray, error **err);
size_t rc_get_real_array_asfloat(confFile *rc, char* key, float** parray, error **err);
size_t rc_get_integer_asint(confFile *rc, char* key, int** parray, error **err);
size_t rc_get_string_array(confFile *rc, char* key, char*** parray, error **err);
size_t rc_safeget_real_array(confFile *rc, char* key, double** parray, error **err);
size_t rc_safeget_integer_array(confFile *rc, char* key, long** parray, error **err);
size_t rc_safeget_real_array_asfloat(confFile *rc, char* key, float** parray, error **err);
size_t rc_safeget_integer_array_asint(confFile *rc, char* key, int** parray, error **err);
size_t rc_safeget_string_array(confFile *rc, char* key, char*** parray, error **err);
size_t rc_get_real_asarray(confFile *rc, char* key, double** parray,int sz, error **err);
size_t rc_get_real_assquarematrix(confFile *rc, char *key, double** res,int sz ,error **err);
size_t rc_get_integer_asarray(confFile *rc, char* key, long** parray,int sz, error **err);
size_t rc_get_string_asarray(confFile *rc, char* key, char*** parray,int sz, error **err);

void rc_do_args(confFile *rc, int argc, char **argv, char shortopt, char* longopt, error **err);
void rc_close(confFile **self);
void rc_toggle_display(confFile *rc, FILE* display, error **err);

int rc_get_list(confFile *rc, int npars, char* prefix, int failflag, error **err, ...);

void rc_free(confFile *rc,void** smthng,error **err);

confFile* rc_alias(confFile *rc, char* prefix, error **err);
void rc_add_prefix(char* name, char* prefix, char* res, error **err);
confFile* rc_get_rootrc(confFile *rc, char** pprefix, error **err);

confFile* rc_init_from_args(int argc, char** argv, error **err);

void rc_add_mem(confFile *rc, void* ptr, error **err);
void* rc_add_malloc(confFile *rc, size_t sz, error **err);

#define testErrorRetLUA(test,lua_error_type,L,message,err,line,return_value,...) {\
if (test) {	\
err=addError(lua_error_type,(char*)lua_tostring(L,-1),err,line); \
err=topErrorVA(forwardErr,message,err,line,__VA_ARGS__); \
return return_value; \
} \
}

#endif


#define rc_base		        -600
#define rc_size_error	    -1 + io_base
#define rc_lua_error 	    -2 + io_base
#define rc_maxdepth_error	-3 + io_base
#define rc_bad_type     	-4 + io_base
#define rc_bad_key      	-5 + io_base



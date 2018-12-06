/*
 *  readConf.c
 *  readConf
 *
 *  Created by Karim Benabed on 15/03/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#include "readConf.h"


confFile* rc_init(char* fnm, error **err) {
  confFile *rc;
  
  rc = malloc_err(sizeof(confFile),err);
  forwardError(*err,__LINE__,NULL);
  
  rc->mem_sz = RC_MEM_STEP;
  rc->mem = malloc_err(sizeof(void*)*rc->mem_sz,err);
  forwardError(*err,__LINE__,NULL);
  rc->mem_cur = 0;
  
  /*rc->bf = malloc_err(sizeof(char)*rc->bf_sz,err);
   forwardError(*err,__LINE__,NULL);
   rc->bf_sz = RC_BFSTEP;
   rc->bf_cur = 0;*/
  
  // save the name of the file
  testErrorRetVA(strlen(fnm)>=RC_FNMLEN,rc_size_error,"Size of filename too long (found %d)",*err,__LINE__,NULL,strlen(fnm));
  sprintf(rc->fnm,"%s",fnm);
  
  rc->last_key[0]='\0';
  rc->display = NULL;
  rc->root[0]='\0';
  
  // initialize the lua interpretor
  rc->L = luaL_newstate();
  testErrorRet(rc->L==NULL,rc_lua_error,"Cannot start LUA",*err,__LINE__,NULL);
  luaL_openlibs(rc->L);
  
  // init special object rc
  lua_newtable(rc->L);
  lua_setglobal(rc->L,"rc");
  
  rc->alias = NULL;
  rc->alias_prefix[0]='\0';
  
  return rc;
}

confFile* rc_open(char* fnm, error **err) {
  confFile *rc;
  int lua_err;
  
  rc = rc_init(fnm,err);
  forwardError(*err,__LINE__,NULL);
  
  // read file
  lua_err = luaL_dofile(rc->L, fnm);
  testErrorRetLUA(lua_err!=0,lua_err,rc->L,"Cannot read file %s",*err,__LINE__,NULL,fnm);
    
  rc->last_key[0]='\0';
  rc->display = NULL;
  
  return rc;
}

confFile* rc_open_with(char* fnm, char* with,int isfile, error **err) {
  confFile *rc;
  int lua_err;
  char _funcname[RC_FNMLEN],_with[RC_FNMLEN];
  char * funcname;
  
  sprintf(_funcname,RC_DEFAULT_FUNCNAME);
  funcname = _funcname;
  
  if ((with==NULL) || (with[0]=='\0')) {
    rc = rc_open(fnm,err);
    forwardError(*err,__LINE__,NULL);
    
  }
  rc = rc_init(fnm,err);
  forwardError(*err,__LINE__,NULL);
  
  // read the parser command
  if (isfile==1) {
    // locate the name of the function to call (if any)
    strcpy(_with,with);
    funcname = rindex(_with,':');
    if (funcname == NULL) {
      funcname = _funcname;
    } else {
      //split the char
      funcname[0]='\0';
      funcname++;
    }
    lua_err = luaL_dofile(rc->L, _with);
    testErrorRetLUA(lua_err!=0,lua_err,rc->L,"Cannot read file %s",*err,__LINE__,NULL,_with);
  } else {
    lua_err = luaL_dostring(rc->L, with);
    testErrorRetLUA(lua_err!=0,lua_err,rc->L,"Cannot understand command '%s'",*err,__LINE__,NULL,with);
  }
  
  // deal with the parameter file
  
  
  lua_getfield(rc->L, LUA_GLOBALSINDEX, funcname);
  testErrorRetVA(lua_isfunction(rc->L, -1)==0,rc_bad_type,
               "Cannot find parser function %s",
               *err,__LINE__,NULL,funcname);
  lua_pushstring(rc->L, fnm);
  lua_err = lua_pcall(rc->L, 1, 1,0);
  testErrorRetLUA(lua_err!=0,lua_err,rc->L,
                  "Cannot apply parser function on file file %s",
                  *err,__LINE__,NULL,fnm);
  
  if (lua_isnil(rc->L,-1)==1) {
    // everything was imported in the global namespace
    rc->root[0]='\0';
  } else {
    // parameters were imported in a specifica namespace
    // I need to name it
    lua_setfield(rc->L,-1,RC_ROOTNAME);
    sprintf(rc->root,RC_ROOTNAME);
  }  
  
  return rc;
}

confFile* rc_alias(confFile *rc, char* prefix, error **err) {
  confFile *rca,*rcr;
  char *root_prefix;
  int ll, llr;
  
  //_DEBUGHERE_("%s",prefix);
  rcr = rc_get_rootrc(rc,&root_prefix,err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("%s",root_prefix);
  
  ll = strlen(prefix);
  llr = strlen(root_prefix);
  testErrorRet(ll+1+llr>RC_FNMLEN,rc_size_error,"string too long",*err,__LINE__,NULL);
  
  rca = malloc_err(sizeof(confFile),err);
  forwardError(*err,__LINE__,NULL);
  
  
  rca->alias = rcr;
  sprintf(rca->alias_prefix,"%s%s",root_prefix,prefix);
  if(ll+llr>1) {
    if (rca->alias_prefix[ll+llr-1]!='.') {
        rca->alias_prefix[ll+llr]='.';
        rca->alias_prefix[ll+1+llr]='\0';
    }
  }  
  rca->display = rc->display;
  //_DEBUGHERE_("","");
  
  return rca;
}

void rc_add_prefix(char* name, char* prefix, char* res, error **err) {
  int ll,llp;
  ll = strlen(name);
  llp = strlen(prefix);
  testErrorRet(ll+llp>RC_FNMLEN,rc_size_error,"string too long",*err,__LINE__,);
  sprintf(res,"%s%s",prefix,name);
}

confFile* rc_get_rootrc(confFile *rc, char** pprefix, error **err) {

  if (pprefix!=NULL) {
    *pprefix = rc->alias_prefix;
  }
    
  if (rc->alias==NULL) {
    testErrorRet(rc->alias_prefix[0]!='\0',rc_bad_type,"Bad rc state",*err,__LINE__,NULL);
    return rc;
  }
  return rc->alias;
}

void rc_toggle_display(confFile *rc, FILE* display, error **err) {
  rc->display = display;    
}

void rc_close(confFile **prc) {
  confFile *rc;
  int imem;
  
  rc = *prc;
  if (rc->alias!=NULL) {
    free(*prc);
    *prc = NULL;
    return;
  }
  
  lua_close(rc->L);
  for (imem=0;imem<rc->mem_cur;imem++) {
    if (rc->mem[imem]!=NULL) {
      free(rc->mem[imem]);
    }
  }
  free(rc);
  *prc = NULL;
  
}

void rc_add_mem(confFile *rc, void* ptr, error **err) {
  confFile *rcr;
  
  rcr = rc_get_rootrc(rc,NULL,err);
  forwardError(*err,__LINE__,);
  
  //testErrorRet(rc->alias!=NULL,rc_bad_type,"Does not work with alias",*err,__LINE__,);
  
  if (rcr->mem_cur==rcr->mem_sz) {
    rcr->mem = resize_err(rcr->mem, rcr->mem_sz*sizeof(void*), (rcr->mem_sz + RC_MEM_STEP)*sizeof(void*), 1, err);
    forwardError(*err,__LINE__,);
    rcr->mem_sz += RC_MEM_STEP;
  }
  rcr->mem[rcr->mem_cur] = ptr;
  rcr->mem_cur++;
}

void* rc_add_malloc(confFile *rc, size_t sz, error **err) {
  void* res;
  
  res = malloc_err(sz,err);
  forwardError(*err,__LINE__,NULL);
  
  rc_add_mem(rc,res,err);
  forwardError(*err,__LINE__,NULL);
  
  return res;
}

int rc_split_key(char* key,size_t pos[RC_MAX_DEPTH*2],error **err) {
  size_t cur;
  size_t depth;
  
  depth = 0;
  cur=0;
  pos[0]=0;
  while(key[cur]!='\0') {
    //fprintf(stderr,"%s -> %d %c\n",key,cur,key[cur]);
    if (key[cur]=='.') {
      testErrorRetVA(depth>=RC_MAX_DEPTH,rc_maxdepth_error,"Reached max depth for key %s",*err,__LINE__,-1,key);
      pos[depth*2+1] = cur;
      depth++;
      pos[depth*2] = cur+1;
    }
    cur++;
  }
  pos[depth*2+1] = cur;
  depth++;
  testErrorRetVA(depth>=RC_MAX_DEPTH,rc_maxdepth_error,
                 "Reached max depth for key %s",
                 *err,__LINE__,-1,key);
  return depth;
}

void rc_get_key_member(char* key, size_t pos[RC_MAX_DEPTH*2],int id, char *cur_key,error **err) {
  size_t cln;
  
  cln = pos[id*2+1]-pos[id*2];
  testErrorRetVA(cln>=RC_FNMLEN,rc_size_error,
                 "key %s is too long (got %d character expected at most %d)",
                 *err,__LINE__,,key,cln,RC_FNMLEN);
  testErrorRetVA(cln==0,rc_size_error,
                 "key %s has a zero length member",
                 *err,__LINE__,,key);
  memcpy(cur_key, key + pos[id*2], cln);
  cur_key[cln]='\0';
  return;
}

int rc_retrieve_key(confFile *rc,char* key,error **err) {
  size_t pos[RC_MAX_DEPTH*2];
  int depth;
  char cur_key[RC_FNMLEN];
  int lua_index,id;
  confFile* rcr;
  char pref_key[RC_FNMLEN+1];
  char *prefix;
  int lk,iid,isArrayElement;
  
  rcr = rc_get_rootrc(rc,&prefix,err);
  forwardError(*err,__LINE__,0);
  
  rc_add_prefix(key,prefix,pref_key,err);
  forwardError(*err,__LINE__,0);
  
  if (strcmp(rcr->last_key, pref_key) == 0) {
    return 1;
  } else if (rcr->last_key[0]!='\0') {
    lua_pop(rcr->L,1);
  }
    
  depth = rc_split_key(pref_key, pos, err);
  forwardError(*err,__LINE__,0);
  
  lua_index = LUA_GLOBALSINDEX;
  if (rcr->root[0]!='\0') {
    lua_getfield(rcr->L, lua_index,cur_key);
    testErrorRet(lua_isnil(rcr->L,-1)==1,rc_lua_error,"cannot retrieve top object",*err,__LINE__,0); 
    lua_index = -1;
  }
  
  for(id=0;id<depth;id++) {
    // get the element
    rc_get_key_member(pref_key, pos, id, cur_key, err);
    forwardError(*err,__LINE__,-1);
    //fprintf(stderr,"key -> %s\n",cur_key);
    //fprintf(stderr,"idx -> %d (%d)\n",lua_index,LUA_GLOBALSINDEX);
    isArrayElement = 0;
    lk = strlen(cur_key);
    if (cur_key[lk-1] == ']') {
      char *vvv;
      int ik;
      //_DEBUGHERE_("found an array element","");
      for(ik=lk-1;cur_key[ik]!='['&&ik>=0;ik--);
      vvv=cur_key + ik+1;
      iid = atoi(vvv);
      //_DEBUGHERE_("%d",iid);
      cur_key[ik]='\0';
      //_DEBUGHERE_("%s",cur_key);
      isArrayElement = 1;
    }
    lua_getfield(rcr->L, lua_index,cur_key);
    if(lua_isnil(rcr->L,-1)==1) {
      break;
    }
    if (isArrayElement == 1) {
      lua_rawgeti(rcr->L, -1, iid);
      lua_remove(rcr->L,-2);      
    }
    if (lua_index != LUA_GLOBALSINDEX) {
      lua_remove(rcr->L,-2);
    }
    lua_index = -1;
  }
  if (id == depth) {
    // success !
    sprintf(rcr->last_key,"%s",pref_key);
    return 1;
  } else {
    // failure
    lua_pop(rcr->L,1);
    rcr->last_key[0]='\0';
    return 0;
  }
}

int rc_has_key(confFile* rc, char* key, error **err) {
  int success;
  success = rc_retrieve_key(rc, key, err);
  forwardError(*err,__LINE__,-1);  
  return success;
}

double rc_get_real(confFile *rc, char *key, error **err) {
  int success;
  double res;
  confFile* rcr;
  char pref_key[RC_FNMLEN+1];
  char *prefix;
  
  rcr = rc_get_rootrc(rc,&prefix,err);
  forwardError(*err,__LINE__,0);
  
  rc_add_prefix(key,prefix,pref_key,err);
  forwardError(*err,__LINE__,0);
  
  success = rc_retrieve_key(rcr, pref_key, err);
  forwardError(*err,__LINE__,-1);
  testErrorRetVA(success==0,rc_bad_key,
                 "Cannot find %s",
                 *err,__LINE__,-1,pref_key);
  testErrorRetVA(lua_isnumber(rcr->L, -1)==0,rc_bad_type,
                 "key %s : Expecting a double got %s",
                 *err,__LINE__,-1,pref_key,lua_tostring(rcr->L,-1));
  res = lua_tonumber(rcr->L,-1);
  if (rc->display!=NULL) {
    fprintf(rc->display,"ReadConf(%s) %s = %g\n",rcr->fnm,pref_key,res);
  }
  return res;
}

double rc_safeget_real(confFile *rc, char *key, double safegard,error **err) {
  int hk;
  double res;
  
  hk = rc_has_key(rc,key,err);
  forwardError(*err,__LINE__,0);
  if (hk==0) {
    return safegard;
  }
  res = rc_get_real(rc,key,err);
  forwardError(*err,__LINE__,0);
  return res;
}

long rc_get_integer(confFile *rc, char *key, error **err) {
  int success;
  long res;
  confFile* rcr;
  char pref_key[RC_FNMLEN+1];
  char *prefix;
  
  rcr = rc_get_rootrc(rc,&prefix,err);
  forwardError(*err,__LINE__,0);
  
  rc_add_prefix(key,prefix,pref_key,err);
  forwardError(*err,__LINE__,0);
  
  success = rc_retrieve_key(rcr, pref_key, err);
  forwardError(*err,__LINE__,-1);
  testErrorRetVA(success==0,rc_bad_key,
                 "Cannot find %s",
                 *err,__LINE__,-1,pref_key);
  testErrorRetVA(lua_isnumber(rcr->L, -1)==0,rc_bad_type,
                 "key %s : Expecting a number got %s",
                 *err,__LINE__,-1,pref_key,lua_tostring(rcr->L,-1));
  res = lua_tointeger(rcr->L,-1);
  if (rc->display!=NULL) {
    fprintf(rc->display,"ReadConf(%s) %s = %ld\n",rcr->fnm,pref_key,res);
  }
  return res;
}

long rc_safeget_integer(confFile *rc, char *key, long safegard,error **err) {
  int hk;
  long res;
  
  hk = rc_has_key(rc,key,err);
  forwardError(*err,__LINE__,0);
  if (hk==0) {
    return safegard;
  }
  res = rc_get_integer(rc,key,err);
  forwardError(*err,__LINE__,0);
  return res;
}


size_t rc_add_string(confFile *rc, char** pbf, size_t *pbf_sz, size_t *pbf_cur, error **err) {
  size_t len,nrsz,bf_sz,bf_cur,delta_res;
  char *res,*bf;

  testErrorRet(rc->alias!=NULL,rc_bad_type,"Does not work with alias",*err,__LINE__,0);
  
  bf_sz = *pbf_sz;
  bf_cur = *pbf_cur;
  bf = *pbf;
  
  len = lua_objlen(rc->L,-1);
  if (bf_sz < bf_cur + len + 1) {
    nrsz = ((bf_cur + len + 1 - bf_sz) / RC_BFSTEP) + 1;
    bf = resize_err(bf, bf_sz , bf_sz + nrsz * RC_BFSTEP,0,err);
    forwardError(*err,__LINE__,-1);
    bf_sz += nrsz * RC_BFSTEP;
  }
  res = bf + bf_cur;
  delta_res = bf_cur;
  memcpy(res, lua_tolstring(rc->L,-1,&len), len);
  res[len] = '\0';
  bf_cur += len + 1;
  
  *pbf_sz = bf_sz;
  *pbf_cur = bf_cur;
  *pbf = bf;
  
  return delta_res;
  
}

char* rc_get_string(confFile *rc, char *key, error **err) {
  int success;
  size_t len;
  char* res;
  confFile* rcr;
  char pref_key[RC_FNMLEN+1];
  char *prefix;
  
  rcr = rc_get_rootrc(rc,&prefix,err);
  forwardError(*err,__LINE__,0);
  
  rc_add_prefix(key,prefix,pref_key,err);
  forwardError(*err,__LINE__,0);
  
  success = rc_retrieve_key(rcr, pref_key, err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(success==0,rc_bad_key,
                 "Cannot find %s",
                 *err,__LINE__,NULL,pref_key);
  testErrorRetVA(lua_isstring(rcr->L, -1)==0,rc_bad_type,
                 "key %s : Expecting a string got %s",
                 *err,__LINE__,NULL,pref_key,lua_tostring(rcr->L,-1));
  
  /*delta_res = rc_add_string(rc, &rc->bf, &rc->bf_sz, &rc->bf_cur, err);
  forwardError(*err,__LINE__,NULL);
  
  return rc->bf + delta_res;*/
  len = lua_objlen(rcr->L,-1);
  res = rc_add_malloc(rcr,sizeof(char) * (len+1),err);
  forwardError(*err,__LINE__,NULL);
  
  memcpy(res, lua_tolstring(rcr->L,-1,&len), len);
  res[len] = '\0';
  if (rc->display!=NULL) {
    fprintf(rc->display,"ReadConf(%s) %s = \"%s\"\n",rcr->fnm,pref_key,res);
  }
  return res;
}

char* rc_safeget_string(confFile *rc, char *key, char* safegard,error **err) {
  int hk;
  char* res;
  
  hk = rc_has_key(rc,key,err);
  forwardError(*err,__LINE__,0);
  if (hk==0) {
    return safegard;
  }
  res = rc_get_string(rc,key,err);
  forwardError(*err,__LINE__,0);
  return res;
}

int rc_is_array(confFile *rc, char* key, error **err) {
  int success;
  confFile* rcr;
  char pref_key[RC_FNMLEN+1];
  char *prefix;
  
  rcr = rc_get_rootrc(rc,&prefix,err);
  forwardError(*err,__LINE__,0);
  
  rc_add_prefix(key,prefix,pref_key,err);
  forwardError(*err,__LINE__,0);
  
  success = rc_retrieve_key(rcr, pref_key, err);
  forwardError(*err,__LINE__,-1);
  testErrorRetVA(success==0,rc_bad_key,
                 "Cannot find '%s'",
                 *err,__LINE__,-1,pref_key);
  
  success = 1;
  if (lua_istable(rcr->L, -1)==0) {
    success = 0;
  }
  return success;
}

size_t rc_array_size(confFile *rc, char* key, error **err) {
  size_t sz;
  int success;
  confFile* rcr;
  char pref_key[RC_FNMLEN+1];
  char *prefix;
  
  rcr = rc_get_rootrc(rc,&prefix,err);
  forwardError(*err,__LINE__,0);
  
  rc_add_prefix(key,prefix,pref_key,err);
  forwardError(*err,__LINE__,0);
  
  success = rc_retrieve_key(rcr, pref_key, err);
  forwardError(*err,__LINE__,-1);
  testErrorRetVA(success==0,rc_bad_key,
                 "Cannot find %s",
                 *err,__LINE__,-1,pref_key);
  testErrorRetVA(lua_istable(rcr->L, -1)==0,rc_bad_type,
               "key %s is not a table",
               *err,__LINE__,-1,pref_key);
  sz = lua_objlen(rcr->L,-1);
  return sz;
}

size_t rc_get_real_array(confFile *rc, char* key, double** parray, error **err) {
  size_t sz,i;
  double *array;
  confFile* rcr;
  char pref_key[RC_FNMLEN+1];
  char *prefix;
  
  rcr = rc_get_rootrc(rc,&prefix,err);
  forwardError(*err,__LINE__,-1);
  
  rc_add_prefix(key,prefix,pref_key,err);
  forwardError(*err,__LINE__,-1);
  
  sz = rc_array_size(rcr, pref_key, err);
  forwardError(*err,__LINE__,-1);
  
  array = rc_add_malloc(rcr,sizeof(double) * sz,err);
  forwardError(*err,__LINE__,-1);
  
  *parray = array;

  for(i=0;i<sz;i++) {
    lua_rawgeti(rcr->L,-1,i+1);
    testErrorVA(lua_isnumber(rcr->L, -1)==0,rc_bad_type,
                   "key %s : Expecting a number got %s",
                   *err,__LINE__,pref_key,lua_tostring(rcr->L,-1));
    if (isError(*err)) {
      lua_pop(rcr->L,1);
      return -1;
    }
    array[i] = lua_tonumber(rcr->L, -1);
    lua_pop(rcr->L,1);
    if (rc->display!=NULL) {
      fprintf(rc->display,"ReadConf(%s) %s[%ld] = %g\n",rcr->fnm,pref_key,i,array[i]);
    }
  }
  
  return sz;
}

size_t rc_safeget_real_array(confFile *rc, char *key, double** safegard,error **err) {
  int hk;
  size_t res;
  
  hk = rc_has_key(rc,key,err);
  forwardError(*err,__LINE__,-1);
  if (hk==0) {
    return 0;
  }
  res = rc_get_real_array(rc,key,safegard,err);
  forwardError(*err,__LINE__,-1);
  return res;
}

size_t rc_get_real_array_asfloat(confFile *rc, char* key, float** parray, error **err) {
  int sz;
  double* lres;
  int ii;
  
  sz = rc_get_real_array(rc,key,&lres,err);
  forwardError(*err,__LINE__,-1);
  
  *parray = rc_add_malloc(rc,sizeof(float)*sz,err);
  forwardError(*err,__LINE__,-1);
  
  for(ii=0;ii<sz;ii++) {
    (*parray)[ii] = lres[ii];
  }
  
  rc_free(rc,&lres,err);
  forwardError(*err,__LINE__,-1);
  
  return sz;
}

size_t rc_safeget_real_array_asfloat(confFile *rc, char* key, float** safegard, error **err) {
  int hk;
  size_t res;
  
  hk = rc_has_key(rc,key,err);
  forwardError(*err,__LINE__,-1);
  if (hk==0) {
    return 0;
  }
  res = rc_get_real_array_asfloat(rc,key,safegard,err);
  forwardError(*err,__LINE__,-1);
  return res;
}

size_t rc_get_real_asarray(confFile *rc, char *key, double** res,int sz,error **err) {
  int flg;
  size_t rsz;
  double scalar;
  int i;
  
  flg = rc_is_array(rc,key,err);
  forwardError(*err,__LINE__,-1);
  
  if (flg==1) {
    rsz = rc_get_real_array(rc,key,res,err);
    forwardError(*err,__LINE__,-1);
    testErrorRetVA(rsz!=sz,rc_size_error,"bad size for array %s (got %d expected %d)",*err,__LINE__,-1,rsz,sz);
    return rsz;
  }
  scalar = rc_get_real(rc,key,err);
  forwardError(*err,__LINE__,-1);
  
  *res = rc_add_malloc(rc,sizeof(double)*sz,err);
  forwardError(*err,__LINE__,-1);
  for(i=0;i<sz;i++) {
    (*res)[i] = scalar;
  }
  return sz;
}

size_t rc_get_real_assquarematrix(confFile *rc, char *key, double** res,int sz ,error **err) {
  size_t rsz;
  double scalar;
  int i;
  double *mres;
  int isar,ndum;
  char scomp[1000];
  double *mn;
  
  rsz = sz*sz;

  
  mres = rc_add_malloc(rc, rsz*sizeof(double), err);
  forwardError(*err,__LINE__,-1);
  memset(mres,0,sizeof(double)*rsz);
  
  isar = rc_is_array(rc,key,err);
  forwardError(*err,__LINE__,-1);

  if (isar==0) {
    scalar = rc_get_real(rc,key,err);
    forwardError(*err,__LINE__,-1);
    for(i=0;i<sz;i++) {
      mres[i*sz+i] = scalar;
    }
    *res = mres;
    return rsz;
  }
  
  sprintf(scomp,"%s[1]",key);
  isar = rc_is_array(rc,scomp,err);
  forwardError(*err,__LINE__,-1);
  if (isar == 1) {
    // matrix
    for(i=0;i<sz;i++) {
      sprintf(scomp,"%s[%d]",key,i+1);
      ndum = rc_get_real_array(rc,scomp,&mn,err);
      forwardError(*err,__LINE__,-1);
      testErrorRetVA(ndum!=sz,rc_size_error,
                     "'%s[%d]' in %s does not have the right number elements (got %d expected %d)",
                     *err,__LINE__,-1,key,i,rc->alias_prefix,ndum,sz);
      memcpy(mres + i*sz,mn,sizeof(double)*sz);
      rc_free(rc,(void**) &mn,err);
      forwardError(*err,__LINE__,-1);
    }    
  } else {
    ndum = rc_get_real_array(rc,key,&mn,err);
    forwardError(*err,__LINE__,-1);
    if (ndum == sz) {
      // that would be a diagonal !
      for(i=0;i<sz;i++) {
        mres[i*sz+i] = mn[i];
      }
    } else if(ndum == rsz) {
      memcpy(mres,mn,sizeof(double)*rsz);
    }
    testErrorRetVA(ndum!=sz && ndum!=rsz ,rc_size_error,
                   "'%s' in %s does not have the right number elements (got %d expected %d or %d)",
                   *err,__LINE__,-1,key,rc->alias_prefix,ndum,sz, sz*sz);
    rc_free(rc,(void**) &mn,err);
    forwardError(*err,__LINE__,-1);
  }
  *res = mres;
  return rsz;
}

size_t rc_get_integer_array(confFile *rc, char* key, long** parray, error **err) {
  size_t sz,i;
  long *array;
  confFile* rcr;
  char pref_key[RC_FNMLEN+1];
  char *prefix;
  
  rcr = rc_get_rootrc(rc,&prefix,err);
  forwardError(*err,__LINE__,-1);
  
  rc_add_prefix(key,prefix,pref_key,err);
  forwardError(*err,__LINE__,-1);
  
  sz = rc_array_size(rcr, pref_key, err);
  forwardError(*err,__LINE__,-1);
  
  array = rc_add_malloc(rcr,sizeof(long) * sz,err);
  forwardError(*err,__LINE__,-1);
  
  *parray = array;
  
  for(i=0;i<sz;i++) {
    lua_rawgeti(rcr->L,-1,i+1);
    testErrorVA(lua_isnumber(rcr->L, -1)==0,rc_bad_type,
                   "key %s : Expecting a number got '%s'",
                   *err,__LINE__,pref_key,lua_tostring(rcr->L,-1));
    if (isError(*err)) {
      lua_pop(rcr->L,1);
      return -1;
    }
    array[i] = lua_tointeger(rcr->L, -1);
    lua_pop(rcr->L,1);
    if (rc->display!=NULL) {
      fprintf(rc->display,"ReadConf(%s) %s[%ld] = %ld\n",rcr->fnm,pref_key,i,array[i]);
    }
  }
  return sz;
}

size_t rc_safeget_integer_array(confFile *rc, char *key, long** safegard,error **err) {
  int hk;
  size_t res;
  
  hk = rc_has_key(rc,key,err);
  forwardError(*err,__LINE__,-1);
  if (hk==0) {
    return 0;
  }
  res = rc_get_integer_array(rc,key,safegard,err);
  forwardError(*err,__LINE__,-1);
  return res;
}

size_t rc_get_integer_array_asint(confFile *rc, char* key, int** parray, error **err) {
  int sz;
  long* lres;
  int ii;
  
  sz = rc_get_integer_array(rc,key,&lres,err);
  forwardError(*err,__LINE__,-1);
  
  *parray = rc_add_malloc(rc,sizeof(int)*sz,err);
  forwardError(*err,__LINE__,-1);
  
  for(ii=0;ii<sz;ii++) {
    (*parray)[ii] = lres[ii];
  }
  
  rc_free(rc,&lres,err);
  forwardError(*err,__LINE__,-1);
  
  return sz;
}

size_t rc_safeget_integer_array_asint(confFile *rc, char* key, int** safegard, error **err) {
  int hk;
  size_t res;
  
  hk = rc_has_key(rc,key,err);
  forwardError(*err,__LINE__,-1);
  if (hk==0) {
    return 0;
  }
  res = rc_get_integer_array_asint(rc,key,safegard,err);
  forwardError(*err,__LINE__,-1);
  return res;
}

  
size_t rc_get_integer_asarray(confFile *rc, char *key, long** res,int sz,error **err) {
  int flg;
  size_t rsz;
  long scalar;
  int i;
  
  flg = rc_is_array(rc,key,err);
  forwardError(*err,__LINE__,-1);
  
  if (flg==1) {
    rsz = rc_get_integer_array(rc,key,res,err);
    forwardError(*err,__LINE__,-1);
    testErrorRetVA(rsz!=sz,rc_size_error,"bad size for array %s (got %d expected %d)",*err,__LINE__,-1,rsz,sz);
    return rsz;
  }
  scalar = rc_get_integer(rc,key,err);
  forwardError(*err,__LINE__,-1);
  
  *res = rc_add_malloc(rc,sizeof(int)*sz,err);
  forwardError(*err,__LINE__,-1);
  for(i=0;i<sz;i++) {
    (*res)[i] = scalar;
  }
  return sz;
}

size_t rc_get_string_array(confFile *rc, char* key, char*** parray, error **err) {
  size_t sz,i;
  char **array;
  char *bf,*b0;
  size_t *darray;
  size_t bf_sz, bf_cur;
  confFile* rcr;
  char pref_key[RC_FNMLEN+1];
  char *prefix;
  
  rcr = rc_get_rootrc(rc,&prefix,err);
  forwardError(*err,__LINE__,0);
  
  rc_add_prefix(key,prefix,pref_key,err);
  forwardError(*err,__LINE__,0);

  sz = rc_array_size(rcr, pref_key, err);
  forwardError(*err,__LINE__,-1);
  
  darray = malloc_err(sizeof(size_t) * sz,err);
  forwardError(*err,__LINE__,-1);
  
  bf_sz = RC_BFSTEP;
  bf_cur = 0;
  bf = malloc_err(sizeof(char) * bf_sz,err);
  forwardError(*err,__LINE__,-1);
    
  for(i=0;i<sz;i++) {
    lua_rawgeti(rcr->L,-1,i+1);
    testErrorRetVA(lua_isstring(rcr->L, -1)==0,rc_bad_type,
                   "key %s : Expecting a number got %s",
                   *err,__LINE__,-1,pref_key,lua_tostring(rcr->L,-1));
    
    darray[i] = rc_add_string(rcr, &bf, &bf_sz, &bf_cur, err);
    forwardError(*err,__LINE__,-1);
    lua_pop(rcr->L,1);
    if (rc->display!=NULL) {
      fprintf(rc->display,"ReadConf(%s) %s[%ld] = \"%s\"\n",rcr->fnm,pref_key,i,bf+darray[i]);
    }    
  }
  
  array = rc_add_malloc(rcr,sizeof(char) * bf_cur + sizeof(char*) * sz,err);
  forwardError(*err,__LINE__,-1);

  b0 = ((char*) array) + sz * sizeof(char**);
  memcpy(b0, bf, bf_cur);
  for(i=0;i<sz;i++) {
    array[i] = b0 + darray[i];
  }
  
  free(bf);
  free(darray);
  
  *parray = array;
  
  return sz;
}

size_t rc_safeget_string_array(confFile *rc, char *key, char*** safegard,error **err) {
  int hk;
  size_t res;
  
  hk = rc_has_key(rc,key,err);
  forwardError(*err,__LINE__,-1);
  if (hk==0) {
    return 0;
  }
  res = rc_get_string_array(rc,key,safegard,err);
  forwardError(*err,__LINE__,-1);
  return res;
}

size_t rc_get_string_asarray(confFile *rc, char *key, char*** res,int sz,error **err) {
  int flg;
  size_t rsz;
  char *scalar,*delta;
  int i;
  int ll;
  
  flg = rc_is_array(rc,key,err);
  forwardError(*err,__LINE__,-1);
  
  if (flg==1) {
    rsz = rc_get_string_array(rc,key,res,err);
    forwardError(*err,__LINE__,-1);
    testErrorRetVA(rsz!=sz,rc_size_error,"bad size for array %s (got %d expected %d)",*err,__LINE__,-1,rsz,sz);
    return rsz;
  }
  scalar = rc_get_string(rc,key,err);
  forwardError(*err,__LINE__,-1);
  
  ll = strlen(scalar);
  
  *res = rc_add_malloc(rc,sizeof(char*)*sz + (ll+1)*sizeof(char),err);
  forwardError(*err,__LINE__,-1);
  
  delta = ((char*) *res)+sz*sizeof(char*);
  sprintf(delta,"%s",scalar);
  for(i=0;i<sz;i++) {
    (*res)[i] = delta;
  }
  return sz;
}

void rc_do_args(confFile *rc, int argc, char **argv, char shortopt, char* longopt, error **err) {
  // I should (or not) use getops here at some point...
  int i;
  int ind,pos;
  size_t ln;
  int lua_err;
  
  testErrorRet(rc->alias!=NULL,rc_bad_type,"Does not work with alias",*err,__LINE__,);
  
  ln = strlen(longopt);
  
  for(i=0;i<argc;i++) {
    // an options ?
    if (argv[i][0] != '-')
      continue;
    // cas short
    if (argv[i][1]==shortopt) {
      if (argv[i][2]=='\0') {
        ind = i+1;
        pos = 0;
        break;
      }
      if (argv[i][2]=='=') {
        ind = i;
        pos = 3;
        break;
      } 
    }
    // cas long -
    if (strncmp(argv[i]+1,longopt,ln)==0) {
      if (argv[i][ln+1] == '\0') {
        ind = i+1;
        pos = 0;
        break;
      }
      if (argv[i][ln+1]=='=') {
        ind = i;
        pos = ln+2;
        break;
      }
    }
    // cas long --
    if (argv[i][1]=='-' && strncmp(argv[i]+2,longopt,ln)==0) {
      if (argv[i][ln+2] == '\0') {
        ind = i+1;
        pos = 0;
        break;
      }
      if (argv[i][ln+2]=='=') {
        ind = i;
        pos = ln+3;
        break;
      }
    }
  }
  
  if (i==argc) {
    // rien trouve
    return;
  }
  
  //fprintf(stderr,"-->%s<--\n",argv[ind]+pos);
  lua_err = luaL_dostring(rc->L, argv[ind]+pos);
  testErrorRetLUA(lua_err!=0,lua_err,rc->L,"Cannot understand command '%s'",*err,__LINE__,,argv[ind]+pos);
  
  return;
}

int rc_get_list(confFile *rc, int npars, char* prefix, int failflag, error **err, ...) {
  int ipar;
  va_list valist;
  int isarray;
  char *par_name,*par_type;
  int np;
  void* par_ptr;
  int* par_nel;
  int ip,lp,pp;
  int *par_flag;
  int dummy,li;
  confFile* rca,*rcr;
  
  rca = rc_alias(rc,prefix,err);
  forwardError(*err,__LINE__,0);

  rcr = rc_get_rootrc(rca,NULL,err);
  forwardError(*err,__LINE__,0);

  // start reading the extra parameters 
  va_start(valist,err);
  
  np=0;
  for(ipar=0;ipar<npars;ipar++) {
    
    // get the name of the parameter
    par_name = va_arg(valist,char*);
    pp = strlen(par_name);
    if (pp>=RC_FNMLEN) {
      *err=addErrorVA(rc_size_error,"Error while reading %dth parameter (%s %s) name too lenghty",*err,__LINE__,ipar,par_name,par_type);
      va_end(valist);
      return -1; 
    }
    
    // get type
    par_type = va_arg(valist,char*);
    // get pointer
    par_ptr = va_arg(valist,void*);
    // get nelem if type is array
    isarray=0;
    if (par_type[1]=='*') {
      isarray=1;
      par_nel = va_arg(valist,int*);
    } else if (par_type[1]!='\0') {
      *err=addErrorVA(rc_bad_type,"Error while reading %dth parameter (%s %s) unknown type",*err,__LINE__,ipar,par_name,par_type); 
      va_end(valist);
      return -1; 
    }
    // get flag if type failflag<0
    if (failflag<0) {
      par_flag = va_arg(valist,int*);
    } else {
      par_flag = &dummy;
    }
    *par_flag = 1;
    switch (par_type[0]) {
      case 's': case 'S':
        // do something with strings
        if (isarray==1) {
          (*par_nel) = rc_get_string_array(rca, par_name, (char***) par_ptr, err);
        } else {
          *((char**) par_ptr) = rc_get_string(rca,par_name,err);           
        }
        li=__LINE__;
        break;
        
      case 'd': case 'D':
        // do something with integers
        if (isarray==1) {
          (*par_nel) = rc_get_integer_array(rca, par_name, (long**) par_ptr, err);
        } else {
          *((long*) par_ptr) = rc_get_integer(rca,par_name,err); 
        }
        li=__LINE__;
        break;
        
      case 'g': case 'G':
        // do something with reals
        if (isarray==1) {
          (*par_nel) = rc_get_real_array(rca, par_name, (double**) par_ptr, err);
        } else {
          *((double*) par_ptr) = rc_get_real(rca,par_name,err); 
        }
        break;
        
      default:
        *err=addErrorVA(rc_bad_type,
                        "Error while reading %dth parameter (%s %s) unknown type",
                        *err,__LINE__,ipar,par_name,par_type); 
        va_end(valist);
        return -1;
    }
    
    if (isError(*err)) {
      if (failflag==0) {
        va_end(valist);
        forwardError(*err,li,-1);  
      } else {
        if (failflag<0) {
          *par_flag = 0;
        }
        if (rc->display!=NULL) {
          forwardErrorNoReturn(*err,li);
          fprintf(rc->display,"ReadConf(%s) %s absent -- Error traceback:\n",rcr->fnm,par_name);
          fprintf(rc->display, "+---------------------------------------\n");
          printError(rc->display, *err);
          fprintf(rc->display, "+---------------------------------------\n");
        }
        purgeError(err);
      }
    } else {
      np++;
    }
  }
  
  rc_close(&rca);
  
  return np;  
}


void rc_free(confFile *rc,void** smthng,error **err) {
  void* what;
  int imem;
  confFile *rcr;
  
  rcr = rc_get_rootrc(rc,NULL,err);
  forwardError(*err,__LINE__,);

  what = *smthng;
  *smthng=NULL;
  if (what==NULL) 
    return;
  for(imem=0;imem<rcr->mem_cur;imem++) {
    if (rcr->mem[imem]==what) {
      free(what);
      rcr->mem[imem]=rcr->mem[rcr->mem_cur-1];
      rcr->mem_cur--;
      *smthng=NULL;
      return;
    }
  }
  free(what);
}

confFile* rc_init_from_args(int argc, char** argv, error **err) {
  confFile *rc;
  int verb;
  
  testErrorRet(argc==1,-1,"not enough parameters on the command line",
    *err,__LINE__,NULL); 
     
  printf("%s: reading parfile %s\n",__func__,argv[argc-1]);

  // initialize parameter parsing with a file
  rc = rc_open(argv[argc-1],err);
  quitOnError(*err,__LINE__,stderr);

  // look for modifiers for the parameter file on the command line, prepended with -e or --extra
  rc_do_args(rc,argc,argv,'e',"extra",err);
  forwardError(*err,__LINE__,NULL);
  
  // test rc.verbose
  verb = rc_safeget_integer(rc,"rc.verbose",0,err);
  forwardError(*err,__LINE__,NULL);
  if (verb > 0) { //if positive log on stdout
    rc_toggle_display(rc, stdout, err);      
    forwardError(*err,__LINE__,NULL);
  }
  if (verb < 0) { //if negative log on stderr
    rc_toggle_display(rc, stderr, err);      
    forwardError(*err,__LINE__,NULL);
  } 
  
  return rc;
}

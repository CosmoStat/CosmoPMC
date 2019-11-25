/*
 *  nhist.c
 *  likely
 *
 *  Created by Karim Benabed on 24/02/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifdef __PLANCK__
#include "HL2_likely/pmclib/nhist.h"
#else
#include "nhist.h"
#endif

/* histogram support */
nd_histogram * init_nd_histogram(size_t ndim, size_t *nbins, double *limits, error **err) {
  nd_histogram * self;
  size_t i,tdim;
  
  self=malloc_err(sizeof(nd_histogram),err);
  forwardError(*err,__LINE__,NULL);
  
  self->ndim=ndim;
  self->total=0;
  self->nsamples=0;
  
  tdim=1;
  for(i=0;i<ndim;i++) {
    tdim*=nbins[i];
  }
  self->tdim=tdim;
  testErrorRetVA(tdim<=0,nh_allocate,"Negative or zero binsize (%d)",*err,__LINE__,NULL,tdim);
  
  self->buf=malloc_err(ndim* (sizeof(size_t) + 3*sizeof(double)) + tdim*(2*sizeof(double)),err);
  forwardError(*err,__LINE__,NULL);
  
  self->nbins  = self->buf;
  self->limits = (void*)(((char*) self->buf)    + sizeof(size_t)*ndim);
  self->stps   = (void*)(((char*) self->limits) + sizeof(double)*ndim*2);
  self->data   = (void*)(((char*) self->stps)   + sizeof(double)*ndim);
  self->var    = (void*)(((char*) self->data)   + sizeof(double)*tdim);
  
  self->volume=1;
  for(i=0;i<ndim;i++) {
    self->nbins[i]=nbins[i];
    self->limits[i*2]=limits[i*2];
    self->limits[i*2+1]=limits[i*2+1];
    self->stps[i]=(limits[i*2+1]-limits[i*2])/nbins[i];
    self->volume*=self->stps[i];
  }
  for(i=0;i<tdim;i++) {
    self->data[i]=0;
    self->var[i]=0;
  };

  self->isLog=0;

  return self;
}

void free_nd_histogram(nd_histogram **self) {
  free((*self)->buf);
  *self=NULL;
}

nd_histogram * init_nd_loghistogram(size_t ndim, size_t *nbins, double *limits,error **err) {
  double *loglims;
  nd_histogram *self;
  size_t i;
  loglims=malloc_err(ndim*2*sizeof(double),err);
  forwardError(*err,__LINE__,NULL);
  
        
  for(i=0;i<ndim;i++) {
    loglims[i*2]=log(limits[i*2]);
    loglims[i*2+1]=log(limits[i*2+1]);
  }
  self = init_nd_histogram(ndim, nbins, loglims, err);
  free(loglims);
  forwardError(*err,__LINE__,NULL);
  self->isLog=1;
  return self;
}

#define eps 1.0e-15
void acc_histogram(size_t ndim, size_t nsamples, double * params,double *weights, size_t *pidx, nd_histogram *histo, error **err) {
  size_t i,isample,ip,rp,nb,mul,pos,ap;
  double vp,wg;
  double dwg,vwg;
  int valid;
  wg=1;

  //reinit
  for(i=0;i<histo->tdim;i++) {
    histo->data[i]*=histo->nsamples;
    histo->var[i]*=histo->nsamples*(histo->nsamples-1);
  }

  histo->nsamples=nsamples+histo->nsamples;
  dwg=histo->nsamples;
  vwg=1./((dwg-1.)*dwg);
  dwg=1./dwg;
  
  //do estimate
  for (isample=0; isample<nsamples;isample++) {
    mul=1;
    pos=0;
    valid=1;
    for (ap=0;ap<histo->ndim;ap++) {
      ip=histo->ndim-ap-1;
      rp=pidx[ip];
      //fprintf(stderr,"%d %d %d | ",ap,ip,rp);
      vp=params[rp+isample*ndim];
      if ((vp<=histo->limits[ip*2]) || (vp>=histo->limits[ip*2+1])) {
        valid=0;
        break;
      }
      nb=(vp-histo->limits[ip*2])/histo->stps[ip];
      //protect biggest sample
      if(nb==histo->nbins[ip])
        nb--;
      
      pos+=nb*mul;
      mul*=histo->nbins[ip];
    }
    if (valid==1) {
      if (weights!=NULL) {
        wg=weights[isample];
        testErrorRetVA(!isfinite(wg), nh_infnan, "%g found at %d", *err, __LINE__,,wg,isample);
      }
      histo->data[pos]+=wg*dwg;
      histo->total+=wg;
    }
  }
  
  //do variance
  for (isample=0; isample<nsamples;isample++) {
    mul=1;
    pos=0;
    valid=1;
    for (ap=0;ap<histo->ndim;ap++) {
      ip=histo->ndim-ap-1;
      rp=pidx[ip];
      //fprintf(stderr,"%d %d %d | ",ap,ip,rp);
      vp=params[rp+isample*ndim];
      if ((vp<=histo->limits[ip*2]) || (vp>=histo->limits[ip*2+1])) {
        valid=0;
        break;
      }
      nb=(vp-histo->limits[ip*2])/histo->stps[ip];
      pos+=nb*mul;
      mul*=histo->nbins[ip];
    }
    if (valid==1) {
      if (weights!=NULL)
        wg=weights[isample];
      histo->var[pos]+=(wg-histo->data[pos])*(wg-histo->data[pos])*vwg;
    }
  }  
}
#undef eps

nd_histogram * make_histogram(size_t ndim, size_t *pidx, size_t pdim, size_t nsamples, double * params,
			      double *weights, error **err) {

  double *lims,*pmean,*pvar,*bins;
  size_t ap,rp,ap2;
  double vp,w,tot;
  size_t *nbins,isample;
  nd_histogram *self;
  
  lims=malloc_err(sizeof(double)*ndim*2,err);
  forwardError(*err,__LINE__,NULL);
  pmean=calloc_err(ndim,sizeof(double),err);
  forwardError(*err,__LINE__,NULL);
  pvar=calloc_err(ndim*ndim,sizeof(double),err);
  forwardError(*err,__LINE__,NULL);
  nbins=calloc_err(ndim,sizeof(size_t),err);
  forwardError(*err,__LINE__,NULL);
  for(ap=0;ap<ndim;ap++) {
    lims[ap*2]=params[pidx[ap]];
    lims[2*ap+1]=lims[ap*2];
  }
  tot=0;
  for(isample=0;isample<nsamples;isample++) {
    w=1;
    if (weights!=NULL) {
      w=weights[isample];
    }
    tot+=w;   
    for(ap=0;ap<ndim;ap++) {
      rp=pidx[ap];
      ap2=ap*2;
      vp=params[rp+isample*pdim];
      if (vp<lims[ap2]) {
        lims[ap2]=vp;
      }
      if (vp>lims[ap2+1]) {
        lims[ap2+1]=vp;
      }
      pmean[ap]+=vp*w;
    }
  }
  for(ap=0;ap<ndim;ap++) {
    pmean[ap]/=tot;
  }
  
  tot=0;
  for(isample=0;isample<nsamples;isample++) {
    w=1;
    if (weights!=NULL) {
      w=weights[isample];
    }
    tot+=w;
    for(ap=0;ap<ndim;ap++) {
      rp=pidx[ap];
      vp=params[rp+isample*pdim]-pmean[ap];
      for (ap2=ap;ap2<ndim;ap2++) {
        pvar[ap*ndim+ap2]+=w*vp*(params[pidx[ap2]+isample*pdim]-pmean[ap2]);
      }
    }
  }  

  for(ap=0;ap<ndim;ap++) {
    for(ap2=ap;ap2<ndim;ap2++) {
      pvar[ap*ndim+ap2]/=tot;
      pvar[ap2*ndim+ap]=pvar[ap*ndim+ap2];
     }
  }

  bins=compute_optimal_bin(ndim,nsamples,pvar,err);
  if (isError(*err)) {
    free(lims);
    free(pmean);
    free(pvar);    
    free(nbins);    
    forwardError(*err,__LINE__,NULL);
  }
  free(pvar);

  /*
  fprintf(stderr, "optimal #bin =");
  for (ap=0; ap<ndim; ap++) 
    fprintf(stderr, " %.4f", bins[ap]);
  fprintf(stderr, "\n");
  */
  
  for(ap=0;ap<ndim;ap++) {
    nbins[ap]=ceil((lims[2*ap+1]-lims[2*ap])/bins[ap]);
     //bplus=ceil((lims[2*ap+1]-pmean[ap])/bins[ap]-.5);
     //lims[2*ap+1]=pmean[ap]+(bplus+.5)*bins[ap];
     //bmoins=ceil((-lims[2*ap]+pmean[ap])/bins[ap]+.5);
     //lims[2*ap]=pmean[ap]-(bmoins-.5)*bins[ap];
     //nbins[ap]=1+bmoins+bplus;    
     
  }

  free(bins);
  free(pmean);
  
  self = init_nd_histogram(ndim, nbins, lims,err);
  if (isError(*err)) {
    free(lims);
    free(nbins);    
    forwardError(*err,__LINE__,NULL);
  }

  free(lims);
  free(nbins);    
  acc_histogram(pdim, nsamples, params, weights, pidx, self, err);
  if (isError(*err)) {
    free_nd_histogram(&self);
    forwardError(*err,__LINE__,NULL);
  }

  return self;
}

void fill_histogram(size_t ndim, size_t nsamples, double * params,double *weights, size_t *pidx, nd_histogram *histo, error **err) {
  size_t isample,ip,rp,nb,mul,pos,ap;
  double vp,wg;
  double dwg,vwg,lvol;
  int valid;
  double *rvol;
  size_t *rpos;
  wg=1;
  
  rvol=malloc_err(sizeof(double)*nsamples,err);
  forwardError(*err,__LINE__,);
  rpos=malloc_err(sizeof(size_t)*nsamples,err);
  forwardError(*err,__LINE__,);
  
  histo->nsamples=nsamples;
  dwg=1;
  vwg=0;
  //dwg=histo->nsamples;
  //vwg=1./((dwg-1.)*dwg);
  //dwg=1./dwg;
  
  
  //do estimate
  for (isample=0; isample<nsamples;isample++) {
    mul=1;
    pos=0;
    valid=1;
    lvol=1;
    
    if (weights!=NULL)
      wg=weights[isample];
    vwg+=wg*wg;
    
    for (ap=0;ap<histo->ndim;ap++) {
      ip=histo->ndim-ap-1;
      rp=pidx[ip];
      //fprintf(stderr,"%d %d %d | ",ap,ip,rp);
      vp=params[rp+isample*ndim];
      if (histo->isLog==1) {
        vp=log(vp);
      }
      if ((vp<histo->limits[ip*2]) || (vp>histo->limits[ip*2+1])) {
        valid=0;
        pos=(size_t) -1;
        break;
      }
      nb=(vp-histo->limits[ip*2])/histo->stps[ip];
      if (nb==histo->nbins[ip])
        nb--;
      
      if (histo->isLog==1) {
        lvol*=exp((nb+1)*histo->stps[ip])-exp((nb)*histo->stps[ip]);
      }
      pos+=nb*mul;
      mul*=histo->nbins[ip];
    }
    if (histo->isLog!=1)
      rvol[isample]=histo->volume;
    else
      rvol[isample]=lvol;
    lvol=rvol[isample];
    rpos[isample]=pos;
    if (valid==1) {
      histo->data[pos]+=wg*dwg/lvol;
      histo->total+=wg*dwg;
    }
  }
  
  //do variance
  for (isample=0; isample<nsamples;isample++) {
    double tmp;
    mul=1;
    pos=0;
    valid=1;
    lvol=rvol[isample];
    pos=rpos[isample];
    if (pos == (size_t) -1)
      continue;
    if (weights!=NULL)
      wg=weights[isample];
    tmp = (dwg/lvol-histo->data[pos]);
    histo->var[pos]+=wg*tmp*tmp*vwg;
  }
  
  free(rpos);
  free(rvol);
}

nd_histogram * build_histogram(size_t ndim, size_t *pidx, size_t pdim, size_t nsamples, double * params, double *weights, int isLog, error **err) {
  double *lims,*pmean,*pvar,*bins;
  size_t ap,rp,ap2;
  double ovp,vp,w,tot;
  size_t *nbins,isample;
  nd_histogram *self;
  
  lims=malloc_err(sizeof(double)*ndim*2,err);
  forwardError(*err,__LINE__,NULL);
  pmean=calloc_err(ndim,sizeof(double),err);
  forwardError(*err,__LINE__,NULL);
  pvar=calloc_err(ndim*ndim,sizeof(double),err);
  forwardError(*err,__LINE__,NULL);
  nbins=calloc_err(ndim,sizeof(size_t),err);
  forwardError(*err,__LINE__,NULL);
  for(ap=0;ap<ndim;ap++) {
    lims[ap*2]=params[pidx[ap]];
    lims[2*ap+1]=lims[ap*2];
  }
  tot=0;
  for(isample=0;isample<nsamples;isample++) {
    w=1;
    if (weights!=NULL) {
      w=weights[isample];
    }
    tot+=w;   
    for(ap=0;ap<ndim;ap++) {
      rp=pidx[ap];
      ap2=ap*2;
      vp=params[rp+isample*pdim];
      if (vp<lims[ap2]) {
        lims[ap2]=vp;
      }
      if (vp>lims[ap2+1]) {
        lims[ap2+1]=vp;
      }
      if (isLog==1)
        vp=log(vp);
      pmean[ap]+=vp*w;
    }
  }
  testErrorRet(tot==0, math_singularValue, "Division by zero", *err, __LINE__, NULL);
  for(ap=0;ap<ndim;ap++) {
    pmean[ap]/=tot;
  }
  
  tot=0;
  for(isample=0;isample<nsamples;isample++) {
    w=1;
    if (weights!=NULL) {
      w=weights[isample];
    }
    tot+=w;
    for(ap=0;ap<ndim;ap++) {
      rp=pidx[ap];
      vp=params[rp+isample*pdim];
      if (isLog==1)
        vp=log(vp);
      vp=vp-pmean[ap];
      for (ap2=ap;ap2<ndim;ap2++) {
        ovp=params[pidx[ap2]+isample*pdim];
        if (isLog==1)
          ovp=log(ovp);
        pvar[ap*ndim+ap2]+=w*vp*(ovp-pmean[ap2]);
      }
    }
  }  

  for(ap=0;ap<ndim;ap++) {
    for(ap2=ap;ap2<ndim;ap2++) {
      pvar[ap*ndim+ap2]/=tot;
      pvar[ap2*ndim+ap]=pvar[ap*ndim+ap2];
     }
  }
  bins=compute_optimal_bin(ndim,nsamples,pvar,err);
  if (isError(*err)) {
    free(lims);
    free(pmean);
    free(pvar);    
    free(nbins);    
    forwardError(*err,__LINE__,NULL);
  }
  free(pvar);

  for(ap=0;ap<ndim;ap++) {
    if (isLog==1)
      nbins[ap]=ceil((log(lims[2*ap+1])-log(lims[2*ap]))/bins[ap]);
    else
      nbins[ap]=ceil((lims[2*ap+1]-lims[2*ap])/bins[ap]);      
    /*
    bplus=ceil((lims[2*ap+1]-pmean[ap])/bins[ap]-.5);
    lims[2*ap+1]=pmean[ap]+(bplus+.5)*bins[ap];
    bmoins=ceil((-lims[2*ap]+pmean[ap])/bins[ap]+.5);
    lims[2*ap]=pmean[ap]-(bmoins-.5)*bins[ap];
    nbins[ap]=1+bmoins+bplus;    
    */
  }

  free(bins);
  free(pmean);
  
  if (isLog==1)
    self = init_nd_loghistogram(ndim, nbins, lims,err);
  else
    self = init_nd_histogram(ndim, nbins, lims,err);
  if (isError(*err)) {
    free(lims);
    free(nbins);    
    forwardError(*err,__LINE__,NULL);
  }

  free(lims);
  free(nbins);    
  fill_histogram(pdim, nsamples, params, weights, pidx, self, err);
  if (isError(*err)) {
    free_nd_histogram(&self);
    forwardError(*err,__LINE__,NULL);
  }

  return self;
}

double* compute_optimal_bin(size_t ndim,size_t nsamples, double *sigma, error **err) {
  double *res,mul;
  size_t i;
  res=malloc_err(sizeof(double)*ndim,err);
  forwardError(*err,__LINE__,NULL);
  mul=3.5*pow(nsamples,-1./(2.+ndim));
  for (i=0;i<ndim;i++) {
    res[i]=mul*sqrt(sigma[i*ndim+i]);
  }
  return res;
}


  



#ifndef __XCORR_H
#define __XCORR_H

#include "timexec.h"
#include "pmctools/errorlist.h"
#include "param.h"
//#include "mathstuff.h"

typedef struct {
  char datname[1024];
  int Nth, Mz;
  double **xi, *theta, **sigxi, *nbar, *signbar;
} xcorr_state;

typedef struct {
  int Nth, Mz;
  double **eta, **etainv, *A, *gamma;
} xcorr_model;

/*
double *vector(long nl, long nh, error **err);
void free_vector(double *v);
double **matrix(long nrl, long nrh, long ncl, long nch, error **err);
void free_matrix(double **m);
double **matrix_multiply(const double **a, const double **b, int N, error **err);
double **matrix_copy(const double **M, int N, error **err);
*/

void matrix_print(const double **a, int n, FILE *F);
double matrix_inverse(double **C, int N, error **err);
void matrix_lubksb(double **a, int n, int *indx, double b[], error **err);
void matrix_ludcmp(double **a, int n, int *indx, double *d, error **err);
void matrix_unity_test(const double **a, const double **ainv, int N, error **err);


functions_wrapper_t *init_functions_Xcorr(error **err);
void read_from_config_Xcorr(void **state, FILE *F, error **err);
void init_Xcorr(common_like *like, error **err);
void read_xcorr(const char *name, xcorr_state *xc, error **err);
xcorr_model *init_xcorr_model(int Nth, int Mz, error **err);
void free_xcorr_model(xcorr_model *model);
xcorr_model *params_to_xcorr_model(const xcorr_state *xcorr, const double *params, error **err);
double likeli_Xcorr(common_like *like, const double *params, error **err);
double w_true(double theta, double A, double gamma, error **err);
double w_photoz(double theta, int i, int j, const xcorr_model *model, const double *nbar, error **err);
double xcorr_chi2(const xcorr_state *xcorr, const xcorr_model *model, error **err);
void print_Xcorr(FILE *F, void *state, error **err);
void out_xcorr_model(const xcorr_model *model, const double *theta, FILE *F, error **err);



#define xc_base         -5000
#define xc_allocate     -1 + xc_base
#define xc_file         -2 + xc_base
#define xc_range        -3 + xc_base
#define xc_inconsistent -4 + xc_base
#define xc_unknown      -5 + xc_base
#define xc_singular     -6 + xc_base
#define xc_io		-7 + xc_base
#define xc_infnan       -8 + xc_base
#define xc_unity        -9 + xc_base


#endif

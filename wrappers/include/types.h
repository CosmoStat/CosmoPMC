#ifndef __TYPES_H
#define __TYPES_H

/* Data types */
typedef enum {Mvdens,
	      MixMvdens,
              Lensing,
              SNIa,
              Nz,
              CMB,
              CMBDistPrior,
              BAO,
              Xcorr,
              Mring,
              Clustering,
              GalCorr,
	      ClusterMass,
              Bias,
//	      newmodule
} data_t;

// Increase Ndata_t when adding new module:
#define Ndata_t 14

/* Data type identifier strings */
#define sdata_t(i) ( \
 i==Mvdens       ? "Mvdens" : \
 i==MixMvdens    ? "MixMvdens" : \
 i==Lensing      ? "Lensing" : \
 i==SNIa         ? "SNIa" : \
 i==Nz           ? "Nz" : \
 i==CMB          ? "CMB" : \
 i==CMBDistPrior ? "CMBDistPrior" : \
 i==BAO          ? "BAO" : \
 i==Xcorr        ? "Xcorr" : \
 i==Mring        ? "Mring" : \
 i==Clustering   ? "Clustering"  : \
 i==GalCorr      ? "GalCorr" : \
 i==ClusterMass  ? "ClusterMass" : \
 i==Bias         ? "Bias" : \
 "" )

// Add name of new module to the above pre-processor macro:
// i==newmodule  ? "newmodule" : \


#endif

#ifndef AMREX_MLLINOP_F_H_
#define AMREX_MLLINOP_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_mllinop_apply_bc (const int* lo, const int* hi,
                                 amrex_real* phi, const int* philo, const int* phihi,
                                 const int* mask, const int* mlo, const int* mhi,
                                 const int* mask0, const int* m0lo, const int* m0hi,
                                 const int* mask1, const int* m1lo, const int* m1hi,
                                 const int* mask2, const int* m2lo, const int* m2hi,
                                 const int* mask3, const int* m3lo, const int* m3hi,
                                 const int cdir, const int bct, const amrex_real bcl,
                                 const amrex_real* bcval, const int* blo, const int* bhi,
                                 const int maxorder, const amrex_real* dxinv,
                                 const int inhomog, const int nc, const int cross, const int op_order);

#ifdef __cplusplus
}
#endif

#endif

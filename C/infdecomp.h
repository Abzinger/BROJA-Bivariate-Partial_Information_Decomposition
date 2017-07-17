// File:    Tartu_Information_Decomposition/infdecomp.h
// Author:  Dirk Oliver Theis, Abdullah Makkeh http://theory.cs.ut.ee/
//
// Lang:    GNU C 11
// Doc:     in-code
// ***************************************************
//
// DESCRIPTION
//
//
//
//
// DOCUMENTATION
//
// Mostly in-code; general things are explained in the readme file.
//


#ifndef __Tartu_Inf_Decomp__infdecomp_h__
#define __Tartu_Inf_Decomp__infdecomp_h__

typedef struct TID *  TID_Ptr_t;
typedef unsigned int Error_Code_t;     // 0 = no error / ok

extern
const char * infdecomp_error_message(Error_Code_T)

extern
TID_Ptr_t * infdecomp_setup(unsigned n_X, unsigned n_Y, unsigned n_Z, const double p[],
                            Error_Code_t *p_err);
// Parameters
// ----------
// Input
// `n_X               size of the X-set: X = {0,...,n_X-1}
// `n_Y               size of the Y-set: Y = {0,...,n_Y-1}
// `n_Z               size of the Z-set: Z = {0,...,n_Z-1}
// 


#endif //ndef __Tartu_Inf_Decomp__infdecomp_h__

//EOF Tartu_Information_Decomposition/infdecomp.h

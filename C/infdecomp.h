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

extern
struct TID_Global_Stuff
{
     double negative_probability;  // init'd to   -1.e-5
     double probability_gtr_1;     // init'd to    1.00001
}                                                  TID_global_stuff;

typedef enum   TID_Error
{
     TIDerr_OK                      = 0,
     TIDerr_FATAL_out_of_mem,
     TIDerr_Input_sizeXYZ,
     TIDerr_Input_p_neg,
     TIDerr_Input_p_gtr1,
     TIDerr_Input_marginal_neg
}                                          TID_Error_t;

typedef struct TID    * TID_Ptr_t;

extern
const char * infdecomp_error_message(TID_Error_t);

extern
TID_Ptr_t infdecomp_setup(unsigned n_X, unsigned n_Y, unsigned n_Z, const double p[],
                          TID_Error_t *p_err);
// Parameters
// ==========
// Input
// -----
// n_X                     Size of the X-set: X = {0,...,n_X-1}.
// n_Y                     Size of the Y-set: Y = {0,...,n_Y-1}.
// n_Z                     Size of the Z-set: Z = {0,...,n_Z-1}.
// p[]                     Probability distribution defining the marginals.
//                         Array of size n_X*n_Y*n_Z,
//                         access using p[ x + n_X*( y + n_Y*z ) ]
//
// Output
// ------
// p_err                   Error code; p_err->code==0 if there was no error.
//
// Return
// ------
// anonymous struct        Pointer to internal problem data; pass to other functions.
//                         NULL if there was an error (cf. p_err)



#endif //ndef __Tartu_Inf_Decomp__infdecomp_h__

//EOF Tartu_Information_Decomposition/infdecomp.h

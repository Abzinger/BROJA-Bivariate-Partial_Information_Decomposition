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

struct TID_Parameters
{
     double negative_probability;  // init'd to   -1.e-5
     double probability_gtr_1;     // init'd to    1.00001
};

#define TID_Parameters_Default = (struct TID_Parameters){
     .negative_probability = -1.e-5,
     .probability_gtr_1    = 1.00001
};


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

// extern
// const char * infdecomp_error_message(TID_Error_t);

extern
TID_Ptr_t infdecomp_setup(unsigned n_X, unsigned n_Y, unsigned n_Z, const double p[], struct TID_Parameters,
                          TID_Error_t *p_err);
/*******************************************************************************
Parameters
==========
Input
-----
n_X                     Size of the X-set: X = {0,...,n_X-1}.
n_Y                     Size of the Y-set: Y = {0,...,n_Y-1}.
n_Z                     Size of the Z-set: Z = {0,...,n_Z-1}.
p[]                     Probability distribution defining the marginals.
                        Array of size n_X*n_Y*n_Z,
                        access using p[ x + n_X*( y + n_Y*z ) ]

Output
------
p_err                   Error code; p_err->code==0 if there was no error.

Return
------
                        Pointer to internal problem data; pass it to other fns.
                        NULL if there was an error (cf. p_err)
*******************************************************************************/

void infdecomp_free_TID(TID_Ptr_t *p_tidptr);
/*******************************************************************************
Frees the memory allocated by a previous call to infdecomp_setup().
*******************************************************************************/



struct TID_Solution_Info
{
     double max_viol_primal_nonneg_ieq;
     double max_viol_primal_marg_eqn;
     double max_viol_dual_ieq;
     double max_viol_dual_eqn;
     double dist_to_boundary;  // $\min_{xyz} \max(0,q_{xyz})$
     double screwness;         // $\min_{xyz} \max(0,q_{xyz})/\max(0,q_{*yz})$,  with $0/0 := 1$
     double delta;             // bound on distance to optimum: f(x) - L(x,dual)
};



extern
struct TID_Solution_Info infdecomp_check_solution(TID_Ptr_t, const double q[],
                                                  TID_Error_T *p_err);
/*******************************************************************************
Solves a Linear Program or Geometric Program to decide whether q is optimal for
the instance stored in the TID pointer.

Parameters
==========
Input
-----
TID_Ptr_t               Pointer to problem data returned by a previous call to
                        infdecomp_setup(), and not yet free'd.

q[]                     Your guess of an optimal solution to the instance.
                        Array of size n_X*n_Y*n_Z,
                        access using q[ x + n_X*( y + n_Y*z ) ],
                        where n_X, n_Y, n_Z are the same as in the call to
                        infdecomp_setup().

Output
------
p_err                   Error code; p_err->code==0 if there was no error.

Return
------
                        Structure containing insights into the solution.
                        Undefined content, if there was an error.
*******************************************************************************/

extern
struct TID_Solution_Info infdecomp_optimize(TID_Ptr_t, double epsilon,
                                            double q[], TID_Error_T *p_err);
/*******************************************************************************
Solves a Linear Program or Geometric Program to decide whether q is optimal for
the instance stored in the TID pointer.

Parameters
==========
Input
-----
TID_Ptr_t               Pointer to problem data returned by a previous call to
                        infdecomp_setup(), and not yet free'd.

epsilon                 Optimality goal.

Output
------
q[]                     Near optimal solution.
                        Array of size n_X*n_Y*n_Z,
                        access using q[ x + n_X*( y + n_Y*z ) ],
                        where n_X, n_Y, n_Z are the same as in the call to
                        infdecomp_setup().
                        Undefined if there was an error.

p_err                   Error code; p_err->code==0 if there was no error.

Return
------
                        Structure containing insights into the solution,
                        including the actual bound on the optimality:
                        you want delta <= epsilon, but that may be optimistic.
                        Undefined content, if there was an error.
*******************************************************************************/

/*******************************************************************************

Note: The delta value returned by infdecomp_optimize() may be different from
      the one returned by infdecomp_check_solution() for the same solution.
      The reason is that the two methods for determining delta differ.

*******************************************************************************/

#endif //ndef __Tartu_Inf_Decomp__infdecomp_h__

//EOF Tartu_Information_Decomposition/infdecomp.h

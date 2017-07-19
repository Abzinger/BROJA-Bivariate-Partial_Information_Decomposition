// File:    Tartu_Information_Decomposition/infdecomp_p.h
// Author:  Dirk Oliver Theis, Abdullah Makkeh http://theory.cs.ut.ee/
//
// Lang:    GNU C 11
// Doc:     in-code
// ***************************************************
//
// DESCRIPTION
//
// Private header file for infdecomp.c
//
//
// DOCUMENTATION
//
//
//

#include "infdecomp.h"
#include <stdbool.h>

static
inline unsigned atXYZ(int x, int y, int z, int n_X, int n_Y) { return x + n_X*( y + n_Y*z ); }

// Marginal Equations
// ------------------

struct XxY { int x,y; };
struct XxZ { int x,z; };
struct YxZ { int y,z; };
struct X   { int x;   };

struct Marginal_Eqns
{
     int n_X, n_Y, n_Z;
     int nnz_XY, nnz_XZ;
     double      * b_xy;  // Array of length nnz_XY; RHS of marginal xy[i] is b_xy[i].
     double      * b_xz;  // Array of length nnz_XZ; RHS of marginal xz[i] is b_xz[i].
     struct XxY  * xy;    // Array of length nnz_XY; lexicographic, w/ entries belonging to consecutive x values consecutive (``x-minor storage'')...
     struct XxZ  * xz;    // Array of length nnz_XZ; ...in other words, the elements belonging to the same z (or y, resp.) are consecutive.
};
#define Marginal_Eqns__BAD ((struct Marginal_Eqns){-1,-1,-1,-2,-2,NULL,NULL,NULL,NULL})

// Constructor & destructor
static
struct Marginal_Eqns  create_Marginal_Eqns(int n_X, int n_Y, int n_Z, const double p[],  TID_Error_t *); // Could give out-of-memory error
static
void                  free_Marginal_Eqns  (struct Marginal_Eqns *);



// Variable Data & Mapping Variables <-> Triples
// ---------------------------------------------
struct VarData          // maps (x,yz) pairs to vector indices
{
     int  n_YZ;         // # size of YxZ
     int *begin_yz;     // Start of the x-entries in a YZ layer (array of length n_YZ +1)
};
#define VarData__BAD ( (struct VarData){-1,NULL} )  // dummy values

struct TripleData
{                                 // The triple for variable i ∈ { vars.begin_yz[j],vars.begin_yz[j+1]-1 } is given by:  ( x[i], yz[j] )
     struct VarData      vars;
     int                 n_vars;  // == sum_{j=1}^{vars.n_YZ} (  vars.begin_yz[j+1]-vars.begin_yz[j]  )
     struct YxZ        * yz;      // Array of size vars.n_YZ.
     struct X          * x;       // Array of size n_vars.
};
#define TripleData__BAD ((struct TripleData){VarData__BAD,-1,NULL,NULL})

static
struct TripleData create_TripleData_and_VarData(struct Marginal_Eqns, TID_Error_t *); // Could give out-of-memory error.
static
void              free_TripleData_and_VarData(struct TripleData *);


// Non-Zero Structure of the Hessian
// ---------------------------------

struct Hessian_NZ
{
     int nnz;     // Hessian entries are:
     int *k;      //   H_{k[i],ell[i]}
     int *ell;    // for i=0,...,nnz-1
};

#define Hessian_NZ__BAD ( (struct Hessian_NZ){-1,NULL,NULL} )  // dummy values

// Constructor & desctructor:
static
struct Hessian_NZ create_Hess_NZ (const struct VarData, TID_Error_t *p_err);  // may return (e.g.) TIDerr_FATAL_out_of_mem
static
void              free_Hess_NZ   (struct Hessian_NZ *);


// Full Problem Data
// -----------------

// This is what the Mosek callback functions need to access.
struct TID
{
     struct Marginal_Eqns  margls;
     struct TripleData     triples;
     struct Hessian_NZ     hess_NZ;
};


// TODO: Constructor, desctructor



// Oracle Access to -H(X|YZ)
// -------------------------

static
long double eval    (const struct VarData, const double q[]);                // Reads VarData & q, returns value.

static
void        grad    (const struct VarData, const double q[], double g[]);    // Reads VarData & q, writes g (must be alloc'd)
                                                                             // WARNING: Gives DIV_BY_ZERO unless q[i]>0 for all i !!  (-> use subgrad)
static
bool        subgrad (const struct VarData, const double q[], double sg[]);   // Reads VarData & q, writes sg (must be alloc'd)
                                                                             // WARNING: if q[xyz]=0 but q[*yz]>0 there's no subgradient, and bullshit is taken.
                                                                             // Return: true, if a subgradient exists for q; false, if bullshit had to be used.
extern double bullshit_subgradient;  // Initialized to  -1.e-123 or whatever.

static
void        hess    (const struct VarData, const double q[], double H[]);    // Reads VarData & q, writes H (must be alloc'd);
                                                                             // follows Hessian_NZ structure.
                                                                             // WARNING: Gives Div-By-0 unless q[i]>0 for all i


//EOF Tartu_Information_Decomposition/infdecomp_p.h

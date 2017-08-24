// File:    Tartu_Information_Decomposition/infdecomp.c
// Author:  Dirk Oliver Theis http://theory.cs.ut.ee/
//
// Lang:    *GNU* C 11  [I mean it!]
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
//
//

#include "infdecomp_p.h"
#include <math.h>
#include <stdlib.h>

// ****************************************************
// M a r g i n a l    E q u a t i o n s    S t u f f
//- - - - - - - - - - - - - - - - - - - - - - - - - - -
struct Marginal_Eqns  create_Marginal_Eqns(const int n_X, const int n_Y, const int n_Z, const double p[],
                                           TID_Error_t *p_err)
{
     __label__ FATAL_ERROR;

     struct Marginal_Eqns me = Marginal_Eqns__BAD;
     me.n_X=n_X; me.n_Y=n_Y; me.n_Z=n_Z;

     // XxY - Marginals
     // ---------------

     // Wasting memory here, if many of the RHSs are 0...
     me.xy   = malloc( n_X*n_Y*sizeof(struct XxY) );   if(!me.xy)   { *p_err=TIDerr_FATAL_out_of_mem; goto FATAL_ERROR; }
     me.b_xy = malloc( n_X*n_Y*sizeof(double) );       if(!me.b_xy) { *p_err=TIDerr_FATAL_out_of_mem; goto FATAL_ERROR; }

     me.nnz_XY = 0;
     for(int y=0; y<n_Y; ++y) {
          for(int x=0; x<n_X; ++x) { // ``x-minor storage''
               long double sum = 0.;
               for(int z=0; z<n_Z; ++z)  sum += p[atXYZ(x,y,z,n_X,n_Y)];
               if (sum > 0.) {
                    me.b_xy[me.nnz_XY] = (double)sum;
                    me.xy[me.nnz_XY]   = (struct XxY){.x=x,.y=y};
                    ++me.nnz_XY;
               }
               else if (sum < TID_global_stuff.negative_probability) { *p_err=TIDerr_Input_marginal_neg; goto FATAL_ERROR; }
          }//^ for x
     }//^ for y

     // Free up unneeded space:
     me.xy   = realloc(me.xy,   n_X*n_Y*sizeof(struct XxY) );   if(!me.xy)   { *p_err=TIDerr_FATAL_out_of_mem; goto FATAL_ERROR; }
     me.b_xy = realloc(me.b_xy, n_X*n_Y*sizeof(double) );       if(!me.b_xy) { *p_err=TIDerr_FATAL_out_of_mem; goto FATAL_ERROR; }

     // XxZ - Marginals
     // ---------------

     // Wasting memory here, if many of the RHSs are 0.... :-/
     me.xz   = malloc( n_X*n_Z*sizeof(struct XxZ) );   if(!me.xz)   { *p_err=TIDerr_FATAL_out_of_mem; goto FATAL_ERROR; }
     me.b_xz = malloc( n_X*n_Z*sizeof(double) );       if(!me.b_xz) { *p_err=TIDerr_FATAL_out_of_mem; goto FATAL_ERROR; }

     me.nnz_XZ = 0;
     for(int z=0; z<n_Z; ++z) {
          for(int x=0; x<n_X; ++x) {  // ``x-minor storage''
               long double sum = 0.;
               for(int y=0; y<n_Y; ++y)  sum += p[atXYZ(x,y,z,n_X,n_Y)];
               if (sum > 0.) {
                    me.b_xz[me.nnz_XZ] = (double)sum;
                    me.xz[me.nnz_XZ]   = (struct XxZ){.x=x,.z=z};
                    ++me.nnz_XZ;
               }
               else if (sum < TID_global_stuff.negative_probability) { *p_err=TIDerr_Input_marginal_neg; goto FATAL_ERROR; }
          }//^ for x
     }//^ for z

     // Free up unneeded space:
     me.xz   = realloc(me.xz,   n_X*n_Z*sizeof(struct XxZ) );   if(!me.xz)   { *p_err=TIDerr_FATAL_out_of_mem; goto FATAL_ERROR; }
     me.b_xz = realloc(me.b_xz, n_X*n_Z*sizeof(double) );       if(!me.b_xz) { *p_err=TIDerr_FATAL_out_of_mem; goto FATAL_ERROR; }

     // Done
     // ----
    *p_err = TIDerr_OK;
     return me;


FATAL_ERROR:
     // Tidy up possibly allocated mem.
     // (Keep in mind that `me` is init'd with NULL for the pointers.)
     free(me.xz);
     free(me.xy);
     free(me.b_xz);
     free(me.b_xy);

     return Marginal_Eqns__BAD;
}//^ create_Marginal_Eqns()


void free_Marginal_Eqns(struct Marginal_Eqns * p_me)
{
     free(p_me->xy);
     free(p_me->xz);
     free(p_me->b_xy);
     free(p_me->b_xz);
     *p_me = Marginal_Eqns__BAD;
}//^ free_Marginal_Eqns()

// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// ******************************************************
// V a r D a t a ,    T r i p l e D a t a    s t u f f
//- - - - - - - - - - - - - - - - - - - - - - - - - - - -
struct count_vars { int xyz, yz; };

static void dummy_callback(struct count_vars c, int x, int y, int z) {}

static struct count_vars iterate_vars(const struct Marginal_Eqns, typeof(dummy_callback));

struct TripleData create_TripleData_and_VarData(const struct Marginal_Eqns me, TID_Error_t *p_err)
{
     __label__ FATAL_ERROR;
     auto void addit_callback(struct count_vars, int,int,int);

     struct TripleData triples = TripleData__BAD;

     const struct count_vars n_triples = iterate_vars(me,dummy_callback);

     triples.n_vars    = n_triples.xyz;
     triples.vars.n_YZ = n_triples.yz;

     triples.vars.begin_yz = malloc(triples.vars.n_YZ*sizeof(int));         if(!triples.vars.begin_yz) goto FATAL_ERROR;
     triples.yz            = malloc(triples.vars.n_YZ*sizeof(struct YxZ));  if(!triples.yz)            goto FATAL_ERROR;
     triples.x             = malloc(triples.n_vars   *sizeof(struct X));    if(!triples.x)             goto FATAL_ERROR;


     // Do the work
     // -----------

     struct YxZ current_yz = {.y=-1,.z=-1};  // used inside callback

     iterate_vars(me,addit_callback);

     void addit_callback(struct count_vars counter, int x, int y, int z)
     {
          triples.x[counter.xyz] = (struct X  ){.x=x};
          triples.yz[counter.yz] = (struct YxZ){.y=y,.z=z};
          if (current_yz.y!=y || current_yz.z!=z) {
               current_yz = (struct YxZ){.y=y,.z=z};
               triples.vars.begin_yz[counter.yz] = counter.xyz;
          }
     }

     // Done
     // ----
     *p_err = TIDerr_OK;
     return triples;


FATAL_ERROR:
     // Tidy up possibly allocated stuff.
     // (Remember, each of the following pointers was init'd to NULL.)
     free(triples.vars.begin_yz);
     free(triples.yz);
     free(triples.x);

     *p_err=TIDerr_FATAL_out_of_mem;
     return TripleData__BAD;
}//^ create_TripleData_and_VarData()

struct count_vars
iterate_vars(const struct Marginal_Eqns me, typeof(dummy_callback) cbfun)
{
     // Safer coding: I'm const'ing the array's in `me`.
     const struct XxY  * const me_xy = me.xy;
     const struct XxZ  * const me_xz = me.xz;

     struct count_vars counter = {.xyz=0, .yz=0};
     for (int y=0, begin_xy=0; y<me.n_Y; ++y) {
          if (y < me_xy[begin_xy].y) continue; // skip this y

          for (int z=0,i_xz=0; z<me.n_Z; ++z) {
               if (z < me_xz[i_xz].z) continue; // skip this z

               int x = -1;
               int i_xy = begin_xy;
               while (me_xy[i_xy].y==y  &&  me_xz[i_xz].z==z) {
                    while (me_xy[i_xy].x < me_xz[i_xz].x  &&  me_xy[i_xy].y==y)  ++i_xy;
                    while (me_xz[i_xz].x < me_xy[i_xy].x  &&  me_xz[i_xz].z==z)  ++i_xz;
                    if (me_xy[i_xy].y==y  &&  me_xz[i_xz].z==z) {
                         x = me_xy[i_xy].x;
                         cbfun(counter,x,y,z);
                         ++ counter.xyz;
                    }
               }
               if (x>=0) ++counter.yz;

               while (me_xz[i_xz++].z == z); // step to next z
          }//^ for z
          while (me_xy[begin_xy++].y == y);  // step to next y
     }//^ for y
     return counter;
}//^ iterate_vars()


void free_TripleData_and_VarData(struct TripleData *p_td)
{
     free( p_td->vars.begin_yz );
     free( p_td->yz );
     free( p_td->x );
     *p_td = TripleData__BAD;
}//^ free_TripleData_and_VarData()

// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// ******************************************
// extern  S E T U P
//- - - - - - - - - - - - - - - - - - - - - -

TID_Ptr_t infdecomp_setup(unsigned n_X, unsigned n_Y, unsigned n_Z, const double p[],
                          TID_Error_t *p_err)
{
     __label__ FATAL_ERROR;

     struct Marginal_Eqns  margls   = Marginal_Eqns__BAD;
     struct TripleData     triples  = TripleData__BAD;
     struct Hessian_NZ     hess_NZ  = Hessian_NZ__BAD;
     struct TID           *p_tid    = NULL;

     // Check input data
     // ----------------
     if (n_X<2 || n_Y<2 || n_Z<2)                  *p_err = TIDerr_Input_sizeXYZ;
     if (n_X*(double)n_Y*(double)n_Z > 134217727)  *p_err = TIDerr_Input_sizeXYZ;
     if(*p_err) goto FATAL_ERROR;

     for (int z=0; z<n_Z; ++z) {
          for (int y=0; y<n_Y; ++y) {
               for (int x=0; x<n_X; ++x) {
                    if ( p[atXYZ(x,y,z,n_X,n_Y)] < TID_global_stuff.negative_probability ) *p_err  = TIDerr_Input_p_neg;
                    if ( p[atXYZ(x,y,z,n_X,n_Y)] > TID_global_stuff.probability_gtr_1 )    *p_err  = TIDerr_Input_p_gtr1;
               }
          }
     }
     if(*p_err) goto FATAL_ERROR;


     margls   = create_Marginal_Eqns(n_X,n_Y,n_Z, p, p_err);   if(*p_err) goto FATAL_ERROR;
     triples  = create_TripleData_and_VarData(margls,p_err);   if(*p_err) goto FATAL_ERROR;
     hess_NZ  = create_Hess_NZ(triples.vars, p_err);           if(*p_err) goto FATAL_ERROR;

     p_tid = malloc(sizeof(struct TID));                       if(!p_tid) goto FATAL_ERROR;
     *p_tid = (struct TID){.margls=margls,.triples=triples,.hess_NZ=hess_NZ};

     *p_err = TIDerr_OK;
     return p_tid;


FATAL_ERROR:
     free(p_tid);
     free_Hess_NZ(&hess_NZ);
     free_TripleData_and_VarData(&triples);
     free_Marginal_Eqns(&margls);
     return NULL;
}//^ infdecomp_setup()



// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// ******************************************
// O r a c l e s
//- - - - - - - - - - - - - - - - - - - - - -

long double eval(const struct VarData vars, const double q[])
{
     long double sum = 0.;

     for(int yz=0; yz<vars.n_YZ; ++yz) {
          const int begin_yz = vars.begin_yz[yz  ];
          const int end_yz   = vars.begin_yz[yz+1];

          register long double q_ast = 0.;
          for (int xyz=begin_yz; xyz<end_yz; ++xyz)   q_ast += q[xyz];

          for (int xyz=begin_yz; xyz<end_yz; ++xyz)   if(q[xyz]>0) sum += q[xyz] * logl(q[xyz]/q_ast);
     }//^ for yz

     return sum;
}//^ eval()

void grad(const struct VarData vars, const double q[], double g[])
{
     for(int yz=0; yz<vars.n_YZ; ++yz) {
          const int begin_yz = vars.begin_yz[yz  ];
          const int end_yz   = vars.begin_yz[yz+1];

          register long double q_ast = 0.;
          for (int xyz=begin_yz; xyz<end_yz; ++xyz)   q_ast += q[xyz];

          for (int xyz=begin_yz; xyz<end_yz; ++xyz)   g[xyz] = (double)( logl(q[xyz]/q_ast) );
     }//^ for yz
}//^ grad()

double bullshit_subgradient = -1.e-123;
bool subgrad(const struct VarData vars, const double q[], double sg[])
{
     bool success = true;
     const double bullshit = -1.e-123;

     for(int yz=0; yz<vars.n_YZ; ++yz) {
          const int begin_yz = vars.begin_yz[yz  ];
          const int end_yz   = vars.begin_yz[yz+1];

          register long double q_ast = 0.;
          for (int xyz=begin_yz; xyz<end_yz; ++xyz)   q_ast += q[xyz];

          if (q_ast<=0) {
               for(int xyz=begin_yz; xyz<end_yz; ++xyz)   sg[xyz] = (double) logl( 1./(end_yz-begin_yz) );
          }else {
               for (int xyz=begin_yz; xyz<end_yz; ++xyz) {
                    if      (q[xyz] >  0)                 sg[xyz] = (double) logl( q[xyz] / q_ast );
                    else if (q[xyz] <= 0)                 sg[xyz] = bullshit, success=false;
               }
          }//^ else
     }//^ for yz
     return success;
}//^ subgrad()


// *******************************
// H e s s i a n   S t u f f
// - - - - - - - - - - - - - - - -

void free_Hess_NZ(struct Hessian_NZ *p_hess_nz)
{
     free(p_hess_nz->k);
     free(p_hess_nz->ell);
     *p_hess_nz = Hessian_NZ__BAD;
}//^ free_Hess_NZ()


struct Hessian_NZ create_Hess_NZ(const struct VarData vars, TID_Error_t *p_err)
{
     __label__ FATAL_ERROR;

     int nnz, *k=NULL, *ell=NULL;

     nnz = 0;
     for(int yz=0; yz<vars.n_YZ; ++yz) {
          const int blocksz  =  vars.begin_yz[yz+1] - vars.begin_yz[yz];
          nnz += (blocksz+1)*blocksz/2;
     }//^ for yz

     k   = malloc(nnz*sizeof(int));   if(!k)   goto FATAL_ERROR;
     ell = malloc(nnz*sizeof(int));   if(!ell) goto FATAL_ERROR;

     int idx = 0; // H[idx], k[idx], ell[idx]
     for(int yz=0; yz<vars.n_YZ; ++yz) {
          const int begin_yz = vars.begin_yz[yz  ];
          const int end_yz   = vars.begin_yz[yz+1];

          // Diagonal: H_{xyz,xyz}
          for (int xyz=begin_yz; xyz<end_yz; ++xyz)                                         k[idx]=ell[idx]=xyz,      idx++;

          // Below diagonal: H_{xyz,u}   for xyz<u
          for(int xyz=begin_yz; xyz<end_yz; ++xyz) for(int uyz=begin_yz; uyz<xyz; ++uyz)    k[idx]=xyz, ell[idx]=uyz, idx++;
     }//^ for yz

     *p_err = TIDerr_OK;
     return (struct Hessian_NZ){.nnz=nnz,.k=k,.ell=ell};

FATAL_ERROR:
     // Tidy up possibly allocated mem
     free(k);
     free(ell);

     *p_err=TIDerr_FATAL_out_of_mem;
     return Hessian_NZ__BAD;
}//^ create_Hess_NZ()

void hess(const struct VarData vars, const double q[], double H[])
{
     int idx = 0; // H[idx]

     for(int yz=0; yz<vars.n_YZ; ++yz) {
          const int begin_yz = vars.begin_yz[yz  ];
          const int end_yz   = vars.begin_yz[yz+1];

          register long double q_ast = 0.;
          for (int xyz=begin_yz; xyz<end_yz; ++xyz)   q_ast += q[xyz];

          // Diagonal: H_{xyz,xyz}
          for (int xyz=begin_yz; xyz<end_yz; ++xyz)                                         H[idx ++]  = (q_ast-q[xyz])/(q_ast*q[xyz]);

          // Below diagonal: H_{xyz,u}   for xyz<uyz
          for(int xyz=begin_yz; xyz<end_yz; ++xyz) for(int uyz=begin_yz; uyz<xyz; ++uyz)    H[idx ++]  = -1/q_ast;
     }//^ for yz
}//^ grad()


// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// ******************************************
// M o s e k    I n t e r f a c e
//- - - - - - - - - - - - - - - - - - - - - -

#include <mosek.h>

static int MSKAPI log__callback();
static int MSKAPI structure__callback();
static int MSKAPI evaluation__callback();

// This is near copy-paste from the Mosek documentation:
MSKint32t MSKAPI log_callback(MSKtask_t            const mosek_task,
                              MSKuserhandle_t      const my_handle,
                              MSKcallbackcodee     const caller,
                              const MSKrealt     * const douinf,
                              const MSKint32t    * const intinf,
                              const MSKint64t    * const lintinf)
{
     struct TID const * p_tid = (struct TID *)mosek_task;

     const char *txt;

     switch (caller) {

     case MSK_CALLBACK_BEGIN_ROOT_CUTGEN:
          txt = "root cut generation is started.";
          break;
     case MSK_CALLBACK_IM_ROOT_CUTGEN:
          txt = "within root cut generation at an intermediate stage.";
          break;
     case MSK_CALLBACK_END_ROOT_CUTGEN:
          txt = "root cut generation is is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_OPTIMIZER:
          txt = "the optimizer is started.";
          break;
     case MSK_CALLBACK_END_OPTIMIZER:
          txt = "the optimizer is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_PRESOLVE:
          txt = "the presolve is started.";
          break;
     case MSK_CALLBACK_UPDATE_PRESOLVE:
          txt = "within the presolve procedure.";
          break;
     case MSK_CALLBACK_IM_PRESOLVE:
          txt = "within the presolve procedure at an intermediate stage.";
          break;
     case MSK_CALLBACK_END_PRESOLVE:
          txt = "the presolve is completed.";
          break;
     case MSK_CALLBACK_BEGIN_INTPNT:
          txt = "the interior-point optimizer is started.";
          break;
     case MSK_CALLBACK_INTPNT:
          txt = "within the interior-point optimizer after the information database has been updated.";
          break;
     case MSK_CALLBACK_IM_INTPNT:
          txt = "at an intermediate stage within the interior-point optimizer where the information database has not been updated.";
          break;
     case MSK_CALLBACK_END_INTPNT:
          txt = "the interior-point optimizer is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_CONIC:
          txt = "the conic optimizer is started.";
          break;
     case MSK_CALLBACK_CONIC:
          txt = "within the conic optimizer after the information database has been updated.";
          break;
     case MSK_CALLBACK_IM_CONIC:
          txt = "at an intermediate stage within the conic optimizer where the information database has not been updated.";
          break;
     case MSK_CALLBACK_END_CONIC:
          txt = "the conic optimizer is terminated.";
          break;
     case MSK_CALLBACK_PRIMAL_SIMPLEX:
          txt = "within the primal simplex optimizer.";
          break;
     case MSK_CALLBACK_DUAL_SIMPLEX:
          txt = "within the dual simplex optimizer.";
          break;
     case MSK_CALLBACK_BEGIN_BI:
          txt = "The basis identification procedure has been started.";
          break;
     case MSK_CALLBACK_IM_BI:
          txt = "within the basis identification procedure at an intermediate point.";
          break;
     case MSK_CALLBACK_END_BI:
          txt = "the basis identification procedure is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_PRIMAL_BI:
          txt = "within the basis identification procedure when the primal phase is started.";
          break;
     case MSK_CALLBACK_IM_PRIMAL_BI:
          txt = "within the basis identification procedure at an intermediate point in the primal phase.";
          break;
     case MSK_CALLBACK_UPDATE_PRIMAL_BI:
          txt = "within the basis identification procedure at an intermediate point in the primal phase.";
          break;
     case MSK_CALLBACK_END_PRIMAL_BI:
          txt = "within the basis identification procedure when the primal phase is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_DUAL_BI:
          txt = "within the basis identification procedure when the dual phase is started.";
          break;
     case MSK_CALLBACK_IM_DUAL_BI:
          txt = "within the basis identification procedure at an intermediate point in the dual phase.";
          break;
     case MSK_CALLBACK_UPDATE_DUAL_BI:
          txt = "within the basis identification procedure at an intermediate point in the dual phase.";
          break;
     case MSK_CALLBACK_END_DUAL_BI:
          txt = "within the basis identification procedure when the dual phase is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_SIMPLEX_BI:
          txt = "within the basis identification procedure when the simplex clean-up phase is started.";
          break;
     case MSK_CALLBACK_IM_SIMPLEX_BI:
          txt = "within the basis identification procedure at an intermediate point in the simplex clean-up phase. The frequency of the call-backs is controlled by the MSK_IPAR_LOG_SIM_FREQ parameter.";
          break;
     case MSK_CALLBACK_BEGIN_PRIMAL_SIMPLEX_BI:
          txt = "within the basis identification procedure when the primal simplex clean-up phase is started.";
          break;
     case MSK_CALLBACK_UPDATE_PRIMAL_SIMPLEX_BI:
          txt = "within the basis identification procedure at an intermediate point in the primal simplex clean-up phase. The frequency of the call-backs is controlled by the MSK_IPAR_LOG_SIM_FREQ parameter.";
          break;
     case MSK_CALLBACK_END_PRIMAL_SIMPLEX_BI:
          txt = "within the basis identification procedure when the primal clean-up phase is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_PRIMAL_DUAL_SIMPLEX_BI:
          txt = "within the basis identification procedure when the primal-dual simplex clean-up phase is started.";
          break;
     case MSK_CALLBACK_UPDATE_PRIMAL_DUAL_SIMPLEX_BI:
          txt = "within the basis identification procedure at an intermediate point in the primal-dual simplex clean-up phase. The frequency of the call-backs is controlled by the MSK_IPAR_LOG_SIM_FREQ parameter.";
          break;
     case MSK_CALLBACK_END_PRIMAL_DUAL_SIMPLEX_BI:
          txt = "within the basis identification procedure when the primal-dual clean-up phase is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_DUAL_SIMPLEX_BI:
          txt = "within the basis identification procedure when the dual simplex clean-up phase is started.";
          break;
     case MSK_CALLBACK_UPDATE_DUAL_SIMPLEX_BI:
          txt = "within the basis identification procedure at an intermediate point in the dual simplex clean-up phase. The frequency of the call-backs is controlled by the MSK_IPAR_LOG_SIM_FREQ parameter.";
          break;
     case MSK_CALLBACK_END_DUAL_SIMPLEX_BI:
          txt = "within the basis identification procedure when the dual clean-up phase is terminated.";
          break;
     case MSK_CALLBACK_END_SIMPLEX_BI:
          txt = "within the basis identification procedure when the simplex clean-up phase is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_MIO:
          txt = "the mixed-integer optimizer is started.";
          break;
     case MSK_CALLBACK_IM_MIO:
          txt = "at an intermediate point in the mixed-integer optimizer.";
          break;
     case MSK_CALLBACK_NEW_INT_MIO:
          txt = "after a new integer solution has been located by the mixed-integer optimizer.";
          break;
     case MSK_CALLBACK_END_MIO:
          txt = "the mixed-integer optimizer is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_SIMPLEX:
          txt = "the simplex optimizer is started.";
          break;
     case MSK_CALLBACK_BEGIN_DUAL_SIMPLEX:
          txt = "the dual simplex optimizer started.";
          break;
     case MSK_CALLBACK_IM_DUAL_SIMPLEX:
          txt = "at an intermediate point in the dual simplex optimizer.";
          break;
     case MSK_CALLBACK_UPDATE_DUAL_SIMPLEX:
          txt = "in the dual simplex optimizer.";
          break;
     case MSK_CALLBACK_END_DUAL_SIMPLEX:
          txt = "the dual simplex optimizer is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_PRIMAL_SIMPLEX:
          txt = "the primal simplex optimizer is started.";
          break;
     case MSK_CALLBACK_IM_PRIMAL_SIMPLEX:
          txt = "at an intermediate point in the primal simplex optimizer.";
          break;
     case MSK_CALLBACK_UPDATE_PRIMAL_SIMPLEX:
          txt = "in the primal simplex optimizer.";
          break;
     case MSK_CALLBACK_END_PRIMAL_SIMPLEX:
          txt = "the primal simplex optimizer is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_PRIMAL_DUAL_SIMPLEX:
          txt = "the primal-dual simplex optimizer is started.";
          break;
     case MSK_CALLBACK_IM_PRIMAL_DUAL_SIMPLEX:
          txt = "at an intermediate point in the primal-dual simplex optimizer.";
          break;
     case MSK_CALLBACK_UPDATE_PRIMAL_DUAL_SIMPLEX:
          txt = "in the primal-dual simplex optimizer.";
          break;
     case MSK_CALLBACK_END_PRIMAL_DUAL_SIMPLEX:
          txt = "the primal-dual simplex optimizer is terminated.";
          break;
     case MSK_CALLBACK_END_SIMPLEX:
          txt = "the simplex optimizer is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_INFEAS_ANA:
          txt = "the infeasibility analyzer is started.";
          break;
     case MSK_CALLBACK_END_INFEAS_ANA:
          txt = "the infeasibility analyzer is terminated.";
          break;
     case MSK_CALLBACK_IM_PRIMAL_SENSIVITY:
          txt = "at an intermediate stage of the primal sensitivity analysis.";
          break;
     case MSK_CALLBACK_IM_DUAL_SENSIVITY:
          txt = "at an intermediate stage of the dual sensitivity analysis.";
          break;
     case MSK_CALLBACK_IM_MIO_INTPNT:
          txt = "at an intermediate point in the mixed-integer optimizer while running the interior-point optimizer.";
          break;
     case MSK_CALLBACK_IM_MIO_PRIMAL_SIMPLEX:
          txt = "at an intermediate point in the mixed-integer optimizer while running the primal simplex optimizer.";
          break;
     case MSK_CALLBACK_IM_MIO_DUAL_SIMPLEX:
          txt = "at an intermediate point in the mixed-integer optimizer while running the dual simplex optimizer.";
          break;
     case MSK_CALLBACK_BEGIN_PRIMAL_SETUP_BI:
          txt = "the primal BI setup is started.";
          break;
     case MSK_CALLBACK_END_PRIMAL_SETUP_BI:
          txt = "the primal BI setup is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_DUAL_SETUP_BI:
          txt = "the dual BI phase is started.";
          break;
     case MSK_CALLBACK_END_DUAL_SETUP_BI:
          txt = "the dual BI phase is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_PRIMAL_SENSITIVITY:
          txt = "Primal sensitivity analysis is started.";
          break;
     case MSK_CALLBACK_END_PRIMAL_SENSITIVITY:
          txt = "Primal sensitivity analysis is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_DUAL_SENSITIVITY:
          txt = "Dual sensitivity analysis is started.";
          break;
     case MSK_CALLBACK_END_DUAL_SENSITIVITY:
          txt = "Dual sensitivity analysis is terminated.";
          break;
     case MSK_CALLBACK_BEGIN_LICENSE_WAIT:
          txt = "Begin waiting for license.";
          break;
     case MSK_CALLBACK_END_LICENSE_WAIT:
          txt = "End waiting for license.";
          break;
     case MSK_CALLBACK_IM_LICENSE_WAIT:
          txt = "MOSEK is waiting for a license.";
          break;
     case MSK_CALLBACK_BEGIN_QCQO_REFORMULATE:
          txt = "Begin QCQO reformulation.";
          break;
     case MSK_CALLBACK_END_QCQO_REFORMULATE:
          txt = "End QCQO reformulation.";
          break;
     case MSK_CALLBACK_IM_QO_REFORMULATE:
          txt = "at an intermediate stage of the conic quadratic reformulation.";
          break;
     case MSK_CALLBACK_BEGIN_TO_CONIC:
          txt = "Begin conic reformulation.";
          break;
     case MSK_CALLBACK_END_TO_CONIC:
          txt = "End conic reformulation.";
          break;
     case MSK_CALLBACK_BEGIN_FULL_CONVEXITY_CHECK:
          txt = "Begin full convexity check.";
          break;
     case MSK_CALLBACK_END_FULL_CONVEXITY_CHECK:
          txt = "End full convexity check.";
          break;
     case MSK_CALLBACK_IM_FULL_CONVEXITY_CHECK:
          txt = "at an intermediate stage of the full convexity check.";
          break;
     case MSK_CALLBACK_BEGIN_PRIMAL_REPAIR:
          txt = "Begin primal feasibility repair.";
          break;
     case MSK_CALLBACK_END_PRIMAL_REPAIR:
          txt = "End primal feasibility repair.";
          break;
     case MSK_CALLBACK_BEGIN_READ:
          txt = "MOSEK has started reading a problem file.";
          break;
     case MSK_CALLBACK_IM_READ:
          txt = "Intermediate stage in reading.";
          break;
     case MSK_CALLBACK_END_READ:
          txt = "MOSEK has finished reading a problem file.";
          break;
     case MSK_CALLBACK_BEGIN_WRITE:
          txt = "MOSEK has started writing a problem file.";
          break;
     case MSK_CALLBACK_END_WRITE:
          txt = "MOSEK has finished writing a problem file.";
          break;
     case MSK_CALLBACK_READ_OPF_SECTION:
          txt = "A chunk of QQ non-zeros has been read from a problem file.";
          break;
     case MSK_CALLBACK_IM_LU:
          txt = "within the LU factorization procedure at an intermediate point.";
          break;
     case MSK_CALLBACK_IM_ORDER:
          txt = "within the matrix ordering procedure at an intermediate point.";
          break;
     case MSK_CALLBACK_IM_SIMPLEX:
          txt = "within the simplex optimizer at an intermediate point.";
          break;
     case MSK_CALLBACK_READ_OPF:
          txt = "OPF reader.";
          break;
     case MSK_CALLBACK_WRITE_OPF:
          txt = "OPF writer.";
          break;
     case MSK_CALLBACK_SOLVING_REMOTE:
          txt = " the task is being solved on a remote server.";
          break;
     default:
          txt = "Unknown caller.";
     }


     ----------------------------------------------------------------------------------------------------------------------------------------------------

          switch ( caller ) {

          case MSK_CALLBACK_BEGIN_INTPNT:
               printf("Starting interior-point optimizer\n");

               break;

          case MSK_CALLBACK_INTPNT:
               printf("Iterations: %-3d  Time: %6.2f(%.2f)  ",
                      intinf[MSK_IINF_INTPNT_ITER],
                      douinf[MSK_DINF_OPTIMIZER_TIME],
                      douinf[MSK_DINF_INTPNT_TIME]);


               printf("Primal obj.: %-18.6e  Dual obj.: %-18.6e\n",
                      douinf[MSK_DINF_INTPNT_PRIMAL_OBJ],
                      douinf[MSK_DINF_INTPNT_DUAL_OBJ]);

               break;

          case MSK_CALLBACK_END_INTPNT:
               printf("Interior-point optimizer finished.\n");

               break;

          case MSK_CALLBACK_BEGIN_PRIMAL_SIMPLEX:
               printf("Primal simplex optimizer started.\n");

               break;

          case MSK_CALLBACK_UPDATE_PRIMAL_SIMPLEX:
               printf("Iterations: %-3d  ",
                      intinf[MSK_IINF_SIM_PRIMAL_ITER]);

               printf("  Elapsed time: %6.2f(%.2f)\n",
                      douinf[MSK_DINF_OPTIMIZER_TIME],
                      douinf[MSK_DINF_SIM_TIME]);

               printf("Obj.: %-18.6e\n",
                      douinf[MSK_DINF_SIM_OBJ]);

               break;

          case MSK_CALLBACK_END_PRIMAL_SIMPLEX:
               printf("Primal simplex optimizer finished.\n");

               break;

          case MSK_CALLBACK_BEGIN_DUAL_SIMPLEX:
               printf("Dual simplex optimizer started.\n");

               break;

          case MSK_CALLBACK_UPDATE_DUAL_SIMPLEX:
               printf("Iterations: %-3d  ",intinf[MSK_IINF_SIM_DUAL_ITER]);

               printf("  Elapsed time: %6.2f(%.2f)\n",
                      douinf[MSK_DINF_OPTIMIZER_TIME],
                      douinf[MSK_DINF_SIM_TIME]);

               printf("Obj.: %-18.6e\n",douinf[MSK_DINF_SIM_OBJ]);

               break;

          case MSK_CALLBACK_END_DUAL_SIMPLEX:
               printf("Dual simplex optimizer finished.\n");

               break;

          case MSK_CALLBACK_BEGIN_BI:
               printf("Basis identification started.\n");

               break;

          case MSK_CALLBACK_END_BI:
               printf("Basis identification finished.\n");

               break;

          default:
               printf("Unknown caller.\n");
          }



     // return 1 to terminate.
     return 0;
}//^ log_callback()





// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#define TID_Solution_Info__BAD ((struct TID_Solution_Info){ \
                    max_viol_primal_nonneg_ieq = -1.,       \
                    max_viol_primal_marg_eqn = -1.,         \
                    max_viol_dual_ieq = -1.,                \
                    max_viol_dual_eqn = -1.,                \
                    dist_to_boundary = -1.,                 \
                    screwness = -1.,                        \
                    delta = -1.                             \
                    } )






// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//EOF Tartu_Information_Decomposition/infdecomp.c

# infdecomp.jl
module InfDecomp

using InfDecomp_Base
using ExpCone
using GradDesc

using MathProgBase

features_list  = [:Grad,:Jac,:JacVec,:Hess]::Vector{Symbol} # not doing :HessVec right now

# ---------------------------------------------
# M a t h P r o g B a s e    I n t e r f a c e
# ---------------------------------------------

import MathProgBase.initialize
import MathProgBase.features_available
import MathProgBase.eval_f
import MathProgBase.eval_g
import MathProgBase.eval_grad_f
import MathProgBase.jac_structure
import MathProgBase.hesslag_structure
import MathProgBase.eval_jac_g
import MathProgBase.eval_jac_prod
import MathProgBase.eval_jac_prod_t
import MathProgBase.eval_hesslag_prod
import MathProgBase.eval_hesslag
import MathProgBase.isobjlinear
import MathProgBase.isobjquadratic
import MathProgBase.isconstrlinear
# import MathProgBase.obj_expr
# import MathProgBase.constr_expr
import MathProgBase.getreducedcosts
import MathProgBase.getconstrduals
import MathProgBase.getsolution
import MathProgBase.status
import MathProgBase.getsolvetime

# ------------
# B a s i c s
# ------------

# initialize()
function initialize(e::My_Eval, requested_features::Vector{Symbol})
    println("here!")
    for feat in requested_features
        if feat ∉ features_list
            error("infdecomp.jl: initialize():\n-   JuliaOpt:MathProgBase is asking for a feature ($feat) that I don't have.-   Maybe use another solver?")
        end
    end
end


# features_available()
features_available(::My_Eval) = features_list

# Properties:
isobjlinear(::My_Eval)           = false
isobjquadratic(::My_Eval)        = false
isconstrlinear(::My_Eval, ::Integer) = true

#global count_hesseval = 0
# ------------------------------------------
# E v a l u a t i o n :   0 t h   o r d e r
# ------------------------------------------
function eval_f{TFloat_2,TFloat}(e::My_Eval, x::Vector{TFloat_2},dummy::TFloat=Float64(0)) :: TFloat
    local condent::Float64
    if e.TmpFloat==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        condent = condEntropy(e,x,BigFloat(0.))
        setprecision(BigFloat,prc)
    else
        condent = condEntropy(e,x,Float64(0))
    end
    e.count_eval_f += 1
    return -condent
    ;
end #^ eval_f()

# eval_g --- eval of constraint into g
function eval_g(e::My_Eval, g::Vector{Float64}, x::Vector{Float64})  :: Void
    g .= reshape( reshape(x,1,e.n)*e.Gt , e.m, ) .- e.rhs
    e.count_eval_g += 1
    return nothing
    ;
end # eval_g()


# ------------------------------------------
# E v a l u a t i o n :   1 s t   o r d e r
# ------------------------------------------

####################################################################################################
# G L O B A L   V A R I A B L E
global__ill_sol   = "What?!???"
ill_sol__copy_sol = "Oh, no!!!"
function ill_sol__do_copy_sol(x) :: Void
    if global__ill_sol == nothing
        global global__ill_sol = x[:]
    elseif typeof(global__ill_sol)==typeof(x)
        global global__ill_sol .= x
    else
        @show "FUCK!!!!!"
    end
    return nothing
end
function  ill_sol__dont_copy_sol(x) :: Void
    return nothing
end

function set_copy_sol_behaviour(doit::Bool)
    if doit
        global global__ill_sol   = nothing
        global ill_sol__copy_sol = ill_sol__do_copy_sol
    else
        global global__ill_sol   = nothing
        global ill_sol__copy_sol = ill_sol__dont_copy_sol
    end
    ;
end
####################################################################################################


# eval_grad_f --- eval gradient of objective function
function eval_grad_f(e::My_Eval, g::Vector{Float64}, x::Vector{Float64}) :: Void
    if e.TmpFloat==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        ∇f(e,g,x,BigFloat(0.))
        setprecision(BigFloat,prc)
    else
        ∇f(e,g,x,Float64(0))
    end
    e.count_eval_grad_f += 1

    # Crap:
    ill_sol__copy_sol(x)

    return nothing
    ;
end # eval_grad_f()


# Constraint Jacobian
# jac_structure() --- zero-nonzero pattern of constraint Jacobian
jac_structure(e::My_Eval) :: Tuple{Vector{Int64},Vector{Int64}}   =  ( e.Gt_L , e.Gt_K )
# Note: Gt is transposed, so K and L are swapped [K: rows of G^T; L: columns of G^T]

# eval_jac_g() --- constraint Jacobian   -> J
function eval_jac_g(e::My_Eval, J::Vector{Float64}, x::Vector{Float64}) :: Void
    J .= e.Gt.nzval
    e.count_eval_jac_g += 1
    return nothing
    ;
end # eval_jac_g()


# eval_jac_prod() --- constraint_Jacobian * w   -> y
function eval_jac_prod(e::My_Eval, y::Vector{Float64}, x::Vector{Float64}, w::Vector{Float64}) :: Void
    y .= reshape( reshape(w,1,e.n)*e.Gt , e.m, )
    return nothing
    ;
end # eval_jac_prod()


# eval_jac_prod_t() --- constraint_Jacobian^T * w  -> y
function eval_jac_prod_t{T<:AbstractFloat}(e::My_Eval, y::Vector{T}, x::Vector{T}, w::Vector{T}) :: Void
    y .= e.Gt*w
    return nothing
end


# ------------------------------------------
# E v a l u a t i o n :   2 n d   o r d e r
# ------------------------------------------

# Lagrangian:
# L(x, (σ,μ) ) = σ f(x) + μ' G(x)
# Since G(.) is linear, it's Hessian is 0 anyway, so it won't bother us.

# hesslag_structure() --- zero-nonzero pattern of Hessian [wrt x] of the Lagrangian
function hesslag_structure(e::My_Eval)  :: Tuple{Vector{Int64},Vector{Int64}}
    K = Vector{Int64}()
    L = Vector{Int64}()
    counter = 0
    for y = 1:e.n_y
        for z = 1:e.n_z
            # Start with the diagonal
            for x = 1:e.n_x
                i = e.varidx[x,y,z]
                if i>0
                    counter += 1
                    push!(K,i)
                    push!(L,i)
                end
            end
            # Now off-diagonal.
            # H is  treated as a symmetric matrix, but:
            # if both (k,l) & (l,k) are present, their values will be added!
            for x = 1:e.n_x
                for u = 1:(x-1)
                    i_x = e.varidx[x,y,z]
                    i_u = e.varidx[u,y,z]
                    if i_x>0 && i_u>0
                        counter += 1
                        push!(K,i_x)
                        push!(L,i_u)
                    end
                end
            end
        end# for z
    end #^ for y
    return (K,L)
    ;
end # hesslag_structure()


function Hess{TFloat}(e::My_Eval, H::Vector{Float64}, p::Vector{Float64}, σ::Float64, dummy::TFloat) :: Void
    counter = 0
    for y = 1:e.n_y
        for z = 1:e.n_z
            # make marginal P(*yz)
            P_yz  ::TFloat = TFloat(0.)
            for x = 1:e.n_x
                i = e.varidx[x,y,z]
                if i>0
                    P_yz += p[i]
                end
            end

            # now: for all pairs x,u we have:
            # if x ≠  u:   -1/P_yz;
            # if x == u:   ( P_yz - P(xyz) )/(  P_yz * P(xyz) )

            # Start with the diagonal
            for x = 1:e.n_x
                i = e.varidx[x,y,z]
                if i>0
                    counter += 1
                    P_xyz = p[i]
                    H[counter] = Float64(   (P_xyz == 0 ) ?  1.e50  : TFloat(σ)*( P_yz - P_xyz )/(  P_yz * P_xyz )  )
                end
            end
            # Now off-diagonal.
            # H is  treated as a symmetric matrix, but:
            # if both (k,l) & (l,k) are present, their values will be added!
            for x = 1:e.n_x
                for u = 1:(x-1)
                    i_x = e.varidx[x,y,z]
                    i_u = e.varidx[u,y,z]
                    if i_x>0 && i_u>0
                        counter += 1
                        H[counter] = -TFloat(σ)/P_yz
                    end
                end
            end
        end#for z
    end# for y
    return nothing
    ;
end # eval_hesslag()



# eval_hesslag() --- Hessian [wrt x] of the Lagrangian
function eval_hesslag(e::My_Eval, H::Vector{Float64}, x::Vector{Float64}, σ::Float64, μ::Vector{Float64}) :: Void
    if e.TmpFloat==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        Hess(e,H,x,σ,BigFloat(0.))
        setprecision(BigFloat,prc)
        e.count_hesseval += 1
    else
        Hess(e,H,x,σ,Float64(0))
        e.count_hesseval += 1
    end
    return nothing
    ;
end # eval_hesslag()


# eval_hesslag() --- ( Hessian [wrt x] of the Lagrangian ) * v
# function eval_hesslag_prod{T<:AbstractFloat}(e::My_Eval, h::Vector{T}, x::Vector{T}, v::Vector{T}, σ::T, μ::Vector{T}) :: Void
#     h .= σ .* Hf(x)*v
#     return nothing
# end


#------------------
# M A R G I N A L S
#------------------

function marginal_X{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: Array{TFloat_2,1}
    if TFloat_2==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
    end
    marg_x::Array{TFloat_2,1}  = zeros(TFloat_2,e.n_x)
    if TFloat_2==BigFloat
        setprecision(BigFloat,prc)
    end
    for x in 1:e.n_x
        for y in 1:e.n_y
            for z in 1:e.n_z
                i = e.varidx[x,y,z]
                if i>0
                    marg_x[x] += p[i]
                end
            end
        end
    end
    return marg_x
    ;
end#^ marginal_X

function marginal_Y{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: Array{TFloat_2,1}
    if TFloat_2==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
    end
    marg_y::Array{TFloat_2,1}  = zeros(TFloat_2,e.n_y)
    if TFloat_2==BigFloat
        setprecision(BigFloat,prc)
    end
    for y in 1:e.n_y
        for x in 1:e.n_x
            for z in 1:e.n_z
                i = e.varidx[x,y,z]
                if i>0
                    marg_y[y] += p[i]
                end
            end
        end
    end
    return marg_y
    ;
end#^ marginal_Y

function marginal_Z{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: Array{TFloat_2,1}
    if TFloat_2==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
    end
    marg_z::Array{TFloat_2,1}  = zeros(TFloat_2,e.n_z)
    if TFloat_2==BigFloat
        setprecision(BigFloat,prc)
    end
    for z in 1:e.n_z
        for x in 1:e.n_x
            for y in 1:e.n_y
                i = e.varidx[x,y,z]
                if i>0
                    marg_z[z] += p[i]
                end
            end
        end
    end
    return marg_z
    ;
end#^ marginal_Z

function marginal_XY{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: Array{TFloat_2,2}
    if TFloat_2==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        setprecision(BigFloat,prc)
    end
    marg_xy::Array{TFloat_2,2}  = zeros(TFloat_2,e.n_x,e.n_y)
    if TFloat_2==BigFloat
        setprecision(BigFloat,prc)
    end
    for x in 1:e.n_x
        for y in 1:e.n_y
            for z in 1:e.n_z
                i = e.varidx[x,y,z]
                if i>0
                    marg_xy[x,y] += p[i]
                end
            end
        end
    end
    return marg_xy
    ;
end#^ marginals_XY

function marginal_XZ{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: Array{TFloat_2,2}
    if TFloat_2==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
    end
    marg_xz::Array{TFloat_2,2}  = zeros(TFloat_2,e.n_x,e.n_z)
    if TFloat_2==BigFloat
        setprecision(BigFloat,prc)
    end
    for x in 1:e.n_x
        for z in 1:e.n_z
            for y in 1:e.n_y
                i = e.varidx[x,y,z]
                if i>0
                    marg_xz[x,z] += p[i]
                end
            end
        end
    end
    return marg_xz
    ;
end#^ marginal_XZ

function marginal_YZ{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: Array{TFloat_2,2}
    if TFloat_2==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
    end
    marg_yz::Array{TFloat_2,2}  = zeros(TFloat_2,e.n_y,e.n_z)
    if TFloat_2==BigFloat
        setprecision(BigFloat,prc)
    end
    for y in 1:e.n_y
        for z in 1:e.n_z
            for x in 1:e.n_x
                i = e.varidx[x,y,z]
                if i>0
                    marg_yz[y,z] += p[i]
                end
            end
        end
    end
    return marg_yz
    ;
end#^ marginal_YZ


#------------------------------------------
# I N F O R M A T I O N - F U N C T I O N S
#------------------------------------------

# entropy_X --- Shannon entropy of X --- H(X)

function entropy_X{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: TFloat_2
    if TFloat_2==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        s::BigFloat = BigFloat(0.)
        setprecision(BigFloat,prc)
    else
        s = Float64(0)
    end
    marg_x = marginal_X(e, p, TFloat)
    for x in 1:e.n_x
        s  +=  (  ( marg_[x] ≤ 0 )  ?   TFloat(0.)   :   marg_x[x]*log(marg_x[x])  )
    end
    return s
    ;
end#^ entropy_X


# I_X_Y --- Mutual information of  X & Y --- I(X;Y)

function I_X_Y{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: TFloat_2
    if TFloat_2==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        s::BigFloat = BigFloat(0.)
        setprecision(BigFloat,prc)
    else
        s = Float64(0)
    end
    marg_x  = marginal_X(e, p, TFloat)
    marg_y  = marginal_Y(e, p, TFloat)
    marg_xy = marginal_XY(e, p, TFloat)
    for x in 1:e.n_x
        for y in 1:e.n_y
            s  +=  (  (marg_xy[x,y]  ≤ 0 || marg_x[x] ≤ 0 || marg_y[y] ≤ 0)  ?   TFloat_2(0.)   :   marg_xy[x,y]*log( marg_xy[x,y] / (marg_y[y] * marg_x[x]) )  )
        end
    end# ^y
    return s/log(2)
    ;
end#^ I_X_Y


# I_X_Y__Z --- Mutual information of X & Y|Z --- I(X;Y|Z)

function I_X_Y__Z{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: TFloat_2
    if TFloat_2==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        s::BigFloat = BigFloat(0.)
        setprecision(BigFloat,prc)
    else
        s = Float64(0)
    end
    marg_z  = marginal_Z(e, p, TFloat)
    marg_xz = marginal_XZ(e, p, TFloat)
    marg_yz = marginal_YZ(e, p, TFloat)
    for x in 1:e.n_x
        for y in 1:e.n_y
            for z in 1:e.n_z
                i = e.varidx[x,y,z]
                if i>0
                    s  +=  (  (marg_yz[y,z]  ≤ 0 || marg_z[z] ≤ 0 || marg_xz[x,z] ≤ 0 || p[i] ≤ 0 )  ?   TFloat_2(0.)   :   p[i]*log( (p[i] * marg_z[z]) / (marg_xz[x,z] * marg_yz[y,z]) )  )
                end
            end
        end
    end
    return s/log(2)
    ;
end#^ I_X_Y__Z


# I_X_Z__Y --- Mutual information of X & Z|Y --- I(X;Z|Y)

function I_X_Z__Y{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: TFloat_2
    if TFloat_2==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        s::BigFloat = BigFloat(0.)
        setprecision(BigFloat,prc)
    else
        s = Float64(0)
    end
    marg_y  = marginal_Y(e, p, TFloat)
    marg_xy = marginal_XY(e, p, TFloat)
    marg_yz = marginal_YZ(e, p, TFloat)
    for x in 1:e.n_x
        for y in 1:e.n_y
            for z in 1:e.n_z
                i = e.varidx[x,y,z]
                if i>0
                    s  +=  (  (marg_yz[y,z]  ≤ 0 || marg_y[y] ≤ 0 || marg_xy[x,y] ≤ 0 || p[i] ≤ 0 )  ?   TFloat_2(0.)   :   p[i]*log( (p[i] * marg_y[y]) / (marg_xy[x,y] * marg_yz[y,z]) )  )
                end
            end
        end
    end
    return s/log(2)
end#^ I_X_Y__Z


# SI --- Shared Information of Y & Z --- SI(Y;Z)

function SI{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: TFloat_2
    if TFloat_2==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        s::BigFloat = BigFloat(0.)
        setprecision(BigFloat,prc)
    else
        s = Float64(0)
    end
    return I_X_Y(e, p, TFloat) - I_X_Y__Z(e, p, TFloat)
    ;
end#^ SI


# ----------------
# U S I N G   I T
# ----------------

numo_vars(e::My_Eval)                   :: Int64            = e.n
numo_constraints(e::My_Eval)            :: Int64            = e.m
constraints_lowerbounds_vec(e::My_Eval) :: Vector{Float64}  = zeros(Float64,e.m)
constraints_upperbounds_vec(e::My_Eval) :: Vector{Float64}  = zeros(Float64,e.m)
vars_lowerbounds_vec(e::My_Eval)        :: Vector{Float64}  = zeros(Float64,e.n)
vars_upperbounds_vec(e::My_Eval)        :: Vector{Float64}  = [Inf  for j=1:e.n]


using Base.Test


function do_it{T1,T2,T3}(pdf::Dict{Tuple{T1,T2,T3},Float64}, solver, tmpFloat::DataType=BigFloat ;  warmstart::Bool=false, model_type=:IPM)
    print("Starting optimization now (",now(),")!\n")

    # global count_hesseval = 0
    const sd,myeval = create_stuff(pdf,tmpFloat)

    stats = nothing

    if model_type ∈ [:ECOS_L,:SCS_L]
        # Conic Programming, large model
        stats = ExpCone.model_L(myeval,solver)
        return sd,myeval,stats

    elseif model_type ∈ [:ECOS_S,:SCS_S]
        # Conic Programming, small model
        stats = ExpCone.model_S(myeval,solver)
        return sd,myeval,stats

    elseif model_type == :My_GradDesc
        stats = my_gradient_descent(myeval;
                                    max_iter       = 1000000,
                                    eps_grad       = 1.e-20,
                                    eps_steplength = 1.e-20,
                                    stepfactor     = .1)
        return sd,myeval,stats

    else # not Conic Program
        const lb = constraints_lowerbounds_vec(myeval)
        const ub = constraints_upperbounds_vec(myeval)
        const l = vars_lowerbounds_vec(myeval)
        const u = vars_upperbounds_vec(myeval)
        const model = MathProgBase.NonlinearModel(solver)
        MathProgBase.loadproblem!(model, myeval.n, myeval.m, l, u, lb, ub, :Min, myeval)

        if warmstart
            MathProgBase.setwarmstart!(model,ones(Float64,myeval.n))
        end

        MathProgBase.optimize!(model)
        # stat = MathProgBase.status(model)

        return sd,myeval,model;
    end

end #^ do_it()


#-----------------------------------------------
# Solution & Statistics
#-----------------------------------------------
type Solution_and_Stats
    instance_name       :: String
    solver              :: Symbol
    var_num             :: Int64
    x_sz                :: Int64
    y_sz                :: Int64
    z_sz                :: Int64
    status              :: Symbol
    obj_val             :: BigFloat
    q_nonneg_viol       :: BigFloat
    q_min_entry         :: BigFloat
    marginals_1         :: BigFloat
    marginals_2         :: BigFloat
    marginals_Inf       :: BigFloat
    mu_nonneg_viol      :: BigFloat
    complementarity_max :: BigFloat
    complementarity_sum :: BigFloat
    # MI_X_YZ             :: BigFloat
    CI                  :: BigFloat
    SI                  :: BigFloat
    UI_Y                :: BigFloat
    UI_Z                :: BigFloat
    num_eval_f          :: Int64
    num_eval_grad_f     :: Int64
    num_eval_g          :: Int64
    num_eval_jac_g      :: Int64
    num_hessevals       :: Int64
    opt_time            :: BigFloat
    # entropy_X           :: BigFloat
end #^ type Solution_and_Stats


function check_feasibility(name, model, myeval, solver) :: Solution_and_Stats
    if solver ∈ [:ECOS_L,:SCS_L,:ECOS_S,:SCS_S]
        # Conic Program
        stats = model
        model = stats.model

        fstat = Solution_and_Stats(name,solver, 0,0,0,0, status(model) ,  big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),0,0,0,0,0,big(0.))
        fstat.x_sz = myeval.n_x
        fstat.y_sz = myeval.n_y
        fstat.z_sz = myeval.n_z

        fstat.var_num = MathProgBase.numvar(model)
        if status(model) ∈ [:Solve_Succeeded,:Optimal,:NearOptimal,:Suboptimal,:KnitroError,:UserLimit,:FeasibleApproximate,:Error]
            # fstat.obj_val = stats.optimum
            fourtuple = information_quantities(myeval,fstat.obj_val,stats.q)
            (fstat.CI, fstat.SI, fstat.UI_Y, fstat.UI_Z) = fourtuple

            fstat.obj_val = eval_f(myeval,stats.q,Float64(0))

            fstat.q_nonneg_viol = max(-minimum(stats.q),0)
            fstat.q_min_entry   = max(minimum(stats.q),0)

            equation = (stats.q'*myeval.Gt)' - myeval.rhs
            fstat.marginals_1   = norm(equation,1)
            fstat.marginals_2   = norm(equation,2)
            fstat.marginals_Inf = norm(equation,Inf)
        end
        fstat.opt_time = stats.time

        return fstat

    elseif solver == :My_GradDesc
        stats = model

        fstat = Solution_and_Stats(name,solver, 0,0,0,0, stats.status ,  big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),0,0,0,0,0,big(0.))
        fstat.x_sz = myeval.n_x
        fstat.y_sz = myeval.n_y
        fstat.z_sz = myeval.n_z

        fstat.var_num = myeval.n
        if stats.status ∈ [:grad0 :step0]
            # fstat.obj_val = stats.optimum
            fourtuple = information_quantities(myeval,fstat.obj_val,stats.q)
            (fstat.CI, fstat.SI, fstat.UI_Y, fstat.UI_Z) = fourtuple

            fstat.obj_val = eval_f(myeval,stats.q,Float64(0))

            fstat.q_nonneg_viol = max(-minimum(stats.q),0)
            fstat.q_min_entry   = max(minimum(stats.q),0)

            equation = (stats.q'*myeval.Gt)' - myeval.rhs
            fstat.marginals_1   = norm(equation,1)
            fstat.marginals_2   = norm(equation,2)
            fstat.marginals_Inf = norm(equation,Inf)
        end
        fstat.opt_time = stats.time

        return fstat

    else
        fstat = Solution_and_Stats(name,solver, 0,0,0,0, status(model) ,  big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),big(0.),0,0,0,0,0,big(0.))
        fstat.x_sz = myeval.n_x
        fstat.y_sz = myeval.n_y
        fstat.z_sz = myeval.n_z

        fstat.var_num = myeval.n
        if status(model)==:Solve_Succeeded || status(model)==:Optimal || status(model)==:NearOptimal || status(model)==:Suboptimal || status(model)==:KnitroError || ( status(model)==:UserLimit && solver !=:Mosek ) || status(model)==:FeasibleApproximate 
            q = Vector{Float64}( getsolution(model) )
            grad = zeros(Float64, myeval.n)
            InfDecomp.∇f(myeval, grad, q, Float64(.0))

            lambda = Vector{Float64}( getconstrduals(model) )

            mu = grad - myeval.Gt*lambda

            fstat.mu_nonneg_viol = max(-minimum(mu),0)

            fstat.complementarity_max = maximum( abs.(mu) .* abs.(q) )
            fstat.complementarity_sum = sum( abs.(mu) .* abs.(q) )
        else
            if global__ill_sol != nothing
                q = Vector{Float64}( global__ill_sol )
            end
            fstat.mu_nonneg_viol = -10
            fstat.complementarity_max = -10
            fstat.complementarity_sum = -10
        end# ^if status

        fstat.obj_val = eval_f(myeval,q,Float64(0))

        fstat.q_nonneg_viol = max(-minimum(q),0)
        fstat.q_min_entry   = max(minimum(q),0)

        equation = (q'*myeval.Gt)' - myeval.rhs
        fstat.marginals_1   = norm(equation,1)
        fstat.marginals_2   = norm(equation,2)
        fstat.marginals_Inf = norm(equation,Inf)
        # fstat.entropy_X   = entropy_X(myeval,q,Float64(0))
        # fstat.MI_X_YZ     = fstat.entropy_X + fstat.obj_val

        fstat.CI, fstat.SI, fstat.UI_Y, fstat.UI_Z =
            information_quantities(myeval,fstat.obj_val,q)

        fstat.num_eval_f        = myeval.count_eval_f
        fstat.num_eval_grad_f   = myeval.count_eval_grad_f
        fstat.num_eval_g        = myeval.count_eval_g
        fstat.num_eval_jac_g    = myeval.count_eval_jac_g
        fstat.num_hessevals = myeval.count_hesseval
        if solver != :Ipopt
            fstat.opt_time = getsolvetime(model)
        end
        #  end #^if status
        return fstat
    end
    ;
end #^ check_feasibility()

function information_quantities(myeval,optimum,q)
    p = Float64[]
    for x = 1:myeval.n_x
        for y = 1:myeval.n_y
            for z = 1:myeval.n_z
                if myeval.marg_xy[x,y] > 0  &&  myeval.marg_xz[x,z] > 0
                    push!(p,myeval.prb_xyz[x,y,z])
                end
            end# z
        end# y
    end# x

    ci = -(eval_f(myeval,p,Float64(0)) + optimum)/log(2)
    si = SI(myeval,q,Float64(0))
    ui_Y = I_X_Y__Z(myeval,q,Float64(0))
    ui_Z = I_X_Z__Y(myeval,q,Float64(0))

    return ci,si, ui_Y, ui_Z
end

end #^ module

; # EOF

# infdecomp.jl
module InfDecomp

using MathProgBase
# http://mathprogbasejl.readthedocs.io/en/latest/nlp.html

type My_Eval <: MathProgBase.AbstractNLPEvaluator
    n_x     :: Int64
    n_y     :: Int64
    n_z     :: Int64

    n     :: Int64         # number of variables
    m     :: Int64         # number of marginal equations

    varidx  :: Array{Int64,3}                   # 0 if variable not present; otherwise idx of var in {1,...,n}
    xyz     :: Vector{Tuple{Int64,Int64,Int64}} # xyz-triple of varidx

    eqidx   :: Dict{ Tuple{String,Int64,Int64},  Int64} # first idx is "xy" or "xz". Index of eqn in {1,...,m}
    mr_eq   :: Vector{ Tuple{String,Int64,Int64} }      # ("xy",x,y) / ("yz",y,z) of an eqn
    prb_xyz :: Array{Float64,3}
    marg_xy :: Array{Float64,2}
    marg_xz :: Array{Float64,2}

    count_eval_f      :: Int64
    count_eval_grad_f :: Int64
    count_eval_g      :: Int64
    count_eval_jac_g  :: Int64
    count_hesseval    :: Int64

    rhs   :: Vector{Float64}                # m-vector
    Gt    :: SparseMatrixCSC{Float64,Int64} # G^T n x m; transpose of constraint matrix
    Gt_K  :: Vector{Int64}                  # findn()-info of G^T: rows
    Gt_L  :: Vector{Int64}                  # findn()-info of G^T: columns

    TmpFloat       :: DataType
    bigfloat_nbits :: Int64
end

function create_My_Eval(q::Array{Float64,3}, tmpFloat::DataType, bigfloat_precision=256)
    const n_x::Int64 = size(q,1);
    const n_y::Int64 = size(q,2);
    const n_z::Int64 = size(q,3);

    # Create marginals
    prb_xyz::Array{Float64,3}  = zeros(n_x,n_y,n_z)
    marg_xy::Array{Float64,2}  = zeros(n_x,n_y)
    marg_xz::Array{Float64,2}  = zeros(n_x,n_z)

    count_eval_f      :: Int64 = 0
    count_eval_grad_f :: Int64 = 0
    count_eval_g      :: Int64 = 0
    count_eval_jac_g  :: Int64 = 0
    count_hesseval    :: Int64 = 0
    for x in 1:n_x
        for y in 1:n_y
            for z in 1:n_z
                marg_xy[x,y] += q[x,y,z]
                marg_xz[x,z] += q[x,y,z]
            end
        end
    end
    # Find the variables
    varidx::Array{Int64,3}                   = zeros(Bool,size(q));
    xyz   ::Vector{Tuple{Int64,Int64,Int64}} = [ (0,0,0) for i in 1:n_x*n_y*n_z ]
    n::Int64 = 0
    for x in 1:n_x
        for y in 1:n_y
            for z in 1:n_z
                if marg_xy[x,y] > 0  &&  marg_xz[x,z] > 0
                    n += 1
                    varidx[x,y,z] = n
                    xyz[n]        = (x,y,z)
                    prb_xyz[x,y,z]= q[x,y,z]
                else
                    varidx[x,y,z] = 0
                end#if
            end
        end
    end

    # Find the equations
    eqidx = Dict{ Tuple{String,Int64,Int64},Int64}() # first idx is "xy" or "xz"
    mr_eq ::Vector{ Tuple{String,Int64,Int64} }   = [ ("",0,0)   for i in 1:n_x*(n_y+n_z) ]
    m::Int64 = 0
    for x in 1:n_x
        for y in 1:n_y
            if marg_xy[x,y] > 0
                m += 1
                eqidx["xy",x,y] = m
                mr_eq[m]        = ("xy",x,y)
            else
                eqidx["xy",x,y] = 0
            end#if
        end
        for z in 1:n_z
            if marg_xz[x,z] > 0
                m += 1
                eqidx["xz",x,z] = m
                mr_eq[m]        = ("xz",x,z)
            else
                eqidx["xz",x,z] = 0
            end#if
        end
    end #for x

    rhs::Vector{Float64} = zeros(m)
    for k in 1:m
        mr = mr_eq[k]
        if mr[1] == "xy"
            rhs[k] = marg_xy[ mr[2], mr[3] ]
        elseif mr[1]=="xz"
            rhs[k] = marg_xz[ mr[2], mr[3] ]
        else
            print("Fuck! Bug!")
            return;
        end
    end #for all marg equations

    denseGt :: Array{Float64,2} = zeros(n,m)
    for l in 1:n
        (x,y,z) = xyz[l]
        for k in 1:m
            mr = mr_eq[k]
            if mr[1] == "xy"
                xy = mr[2:3]
                if xy[1]==x && xy[2]==y
                    denseGt[l,k] = 1.
                end
            elseif mr[1]=="xz"
                xz = mr[2:3]
                if xz[1]==x && xz[2]==z
                    denseGt[l,k] = 1.
                end
            else
                print("Fuck! Bug!")
                return;
            end
        end
    end

    Gt::SparseMatrixCSC{Float64,Int64} = sparse(denseGt)
    local Gt_K::Array{Int64,1}
    local Gt_L::Array{Int64,1}
    (Gt_K,Gt_L) = findn(Gt)

    TmpFloat       :: DataType  = tmpFloat
    bigfloat_nbits :: Int64     = bigfloat_precision


    return My_Eval(n_x,n_y,n_z, n,m, varidx,xyz, eqidx,mr_eq, prb_xyz, marg_xy,marg_xz, count_eval_f, count_eval_grad_f, count_eval_g, count_eval_jac_g, count_hesseval, rhs, Gt,Gt_K,Gt_L,  TmpFloat,bigfloat_nbits)
    ;
end #^ create_My_eval()


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
function condEntropy{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: TFloat
    # m_yz = marg_yz(e,p,dummy)
    s::TFloat = TFloat(0)
    for y = 1:e.n_y
        for z = 1:e.n_z
            P_yz = TFloat(0)
            for x = 1:e.n_x
                i = e.varidx[x,y,z]
                if i>0
                    P_yz += p[i]
                end
            end
            for x = 1:e.n_x
                i = e.varidx[x,y,z]
                if i>0
                    p_xyz = p[i]
                    s  +=  (  (p_xyz ≤ 0 || P_yz  ≤ 0)  ?   TFloat(0)   :   - p_xyz*log( p_xyz / P_yz )   )
                end
            end
        end
    end
    return s
    ;
end

function eval_f{TFloat_2,TFloat}(e::My_Eval, x::Vector{TFloat_2},dummy::TFloat=Float64(0)) :: TFloat
    local condent::Float64
    if e.TmpFloat==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        condent = condEntropy(e,x,BigFloat(0))
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

function ∇f{TFloat,TFloat_2}(e::My_Eval, grad::Vector{TFloat_2}, p::Vector{TFloat_2}, dummy::TFloat) :: Void
    for y = 1:e.n_y
        for z = 1:e.n_z
            # make marginal P(*yz)
            P_yz::TFloat = TFloat(0.)
            for x = 1:e.n_x
                i = e.varidx[x,y,z]
                if i>0
                    P_yz += p[i]
                end
            end
            # make log-expressions  log( P(xyz) / P(*yz) )
            for x = 1:e.n_x
                i = e.varidx[x,y,z]
                if i>0
                    P_xyz::TFloat = TFloat( p[i] )
                    grad[i] = TFloat_2(   (P_xyz ≤ 0 || P_yz ≤ 0) ?  -log(TFloat(e.n_x))  : log( P_xyz / P_yz )  )
                end
            end
        end# for y
    end# for x
    ;
end # ∇f()

# G L O B A L   V A R I A B L E
ill_sol = nothing
# eval_grad_f --- eval gradient of objective function
function eval_grad_f(e::My_Eval, g::Vector{Float64}, x::Vector{Float64}) :: Void
    if e.TmpFloat==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        ∇f(e,g,x,BigFloat(0))
        setprecision(BigFloat,prc)
    else
        ∇f(e,g,x,Float64(0))
    end
    e.count_eval_grad_f += 1
    # Useful when Mosek status is unknown
    global ill_sol = x      # B U G           <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
        Hess(e,H,x,σ,BigFloat(0))
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


type Set_Data{T1,T2,T3}
    Xidx :: Dict{T1,Int64}
    Yidx :: Dict{T2,Int64}
    Zidx :: Dict{T3,Int64}
    X    :: Vector{T1}
    Y    :: Vector{T2}
    Z    :: Vector{T3}
end

function create_Set_Data{T1,T2,T3}(pdf::Dict{Tuple{T1,T2,T3},Float64}) :: Set_Data{T1,T2,T3}
    Xidx = Dict{T1,Int64}()
    Yidx = Dict{T2,Int64}()
    Zidx = Dict{T3,Int64}()

    anX=aY=aZ=nothing
    for (xyz,val) in pdf
        x,y,z = xyz
        if x ∉ keys(Xidx)
            Xidx[x] = length(Xidx)+1
            anX = x
        end
        if y ∉ keys(Yidx)
            Yidx[y] = length(Yidx)+1
            aY = y
        end
        if z ∉ keys(Zidx)
            Zidx[z] = length(Zidx)+1
            aZ = z
        end
    end

    X ::Vector{T1} = [ anX  for i = 1:length(Xidx) ]
    Y ::Vector{T2} = [ aY   for i = 1:length(Yidx) ]
    Z ::Vector{T3} = [ aZ   for i = 1:length(Zidx) ]

    for (x,i) in Xidx
        X[i] = x
    end
    for (y,i) in Yidx
        Y[i] = y
    end
    for (z,i) in Zidx
        Z[i] = z
    end

    @assert ( length(X) ≥ 2 ) "|Range(X)| ≥ 2 needed"
    @assert ( length(Y) ≥ 2 ) "|Range(Y)| ≥ 2 needed"
    @assert ( length(Z) ≥ 2 ) "|Range(Z)| ≥ 2 needed"

    return Set_Data{T1,T2,T3}(Xidx,Yidx,Zidx,X,Y,Z)
    ;
end #^ create_Set_Data()


function create_stuff{T1,T2,T3}(pdf::Dict{Tuple{T1,T2,T3},Float64}, tmpFloat::DataType) :: Tuple{ Set_Data{T1,T2,T3}, My_Eval }
    sd = create_Set_Data(pdf)

    Q = zeros(Float64,length(sd.X),length(sd.Y),length(sd.Z))
    for (xyz,val) in pdf
        Q[ sd.Xidx[xyz[1]], sd.Yidx[xyz[2]], sd.Zidx[xyz[3]] ] = val
    end

    e = create_My_Eval(Q,tmpFloat)

    return (sd,e)
    ;
end #^ create_stuff

using Base.Test


function do_it{T1,T2,T3}(pdf::Dict{Tuple{T1,T2,T3},Float64}, solver, tmpFloat::DataType=BigFloat, warmstart::Bool=false)
    print("Starting optimization now (",now(),")!\n")

    # global count_hesseval = 0
    const model = MathProgBase.NonlinearModel(solver)
    const sd,myeval = create_stuff(pdf,tmpFloat)
    const lb = constraints_lowerbounds_vec(myeval)
    const ub = constraints_upperbounds_vec(myeval)
    const l = vars_lowerbounds_vec(myeval)
    const u = vars_upperbounds_vec(myeval)

    MathProgBase.loadproblem!(model, myeval.n, myeval.m, l, u, lb, ub, :Min, myeval)

    if warmstart
        MathProgBase.setwarmstart!(model,ones(Float64,myeval.n))
    end

    MathProgBase.optimize!(model)
    stat = MathProgBase.status(model)
    # q = Vector{Float64}( ill_sol )
    # # Check the Reason for Ill_Sol
    # valid = sum(q)
    # max = maximum(q)
    # println("dist ok!, ", valid)
    # for y = 1:myeval.n_y
    #     for z = 1:myeval.n_z
    #         P_yz::tmpFloat = 0
    #         ratio::tmpFloat = 0
    #         i = 0
    #         for x = 1:myeval.n_x
    #             i = myeval.varidx[x,y,z]
    #             if i>0
    #                 P_yz += q[i]
    #                 ratio = P_yz/q[i]
    #                 #               if ( q[i] == sol_boundary && ratio > 1.e10)
    #                 if (ratio > 1.e8)
    #                     println("boundary point found!")
    #                     println("prob =", q[i])
    #                     println("margprob =", P_yz)
    #                     println("scale = ", 1/(max*ratio))
    #                 end
    #             end
    #         end
    #     end# ^z
    # end# ^y
    return sd,myeval,model;

    # println("Status is: ", stat)
    # println("Constraint duals are: ", MathProgBase.getconstrduals(model))
    # println("Reduced costs are: ",  MathProgBase.getreducedcosts(model))
    # println("Solution is: ", MathProgBase.getsolution(model))

    # sub::Array{Int32,1} = [ i for i in 1:myeval.m ]
    # Mosek.getdviolcon(model,1,sub)

    # @test stat == :Optimal
    # @show MathProgBase.getobjval(model)
    # objvaldist = abs( log(n) - MathProgBase.getobjval(model) )*1.e-9
    # println("Distance from true optimum (in 1.e-9): $objvaldist")

    # x = MathProgBase.getsolution(model)
    # dist::Float64 = 0.
    # for j=1:n
    #      dist += abs2(x[j]-1./n)
    # end
    # dist = sqrt(dist)*1.e-9
    # println("Norm Distance from true optimal value (in 1.e-9): $dist")

    # # @test_approx_eq_eps MathProgBase.getobjval(model) log(n)  1.e-300

#    # Test that a second call to optimize! works
#    MathProgBase.setwarmstart!(m,[1,5,5,1])
#    MathProgBase.optimize!(m)
#    stat = MathProgBase.status(m)
#    @test stat == :Optimal
    ;
end #^ do_it()


#-----------------------------------------------
# Solution & Statistics
#-----------------------------------------------
type Solution_and_Stats
    var_num             :: Int64
    x_sz                :: Int64
    y_sz                :: Int64
    z_sz                :: Int64
    status              :: Symbol
    obj_val             :: BigFloat
    q_nonnegativity     :: BigFloat
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


function check_feasibility(model, myeval, solver) :: Solution_and_Stats

    fstat  = Solution_and_Stats( 0,0,0,0, status(model) ,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    fstat.x_sz = myeval.n_x
    fstat.y_sz = myeval.n_y
    fstat.z_sz = myeval.n_z
    fstat.var_num = myeval.n
    if status(model)==:Solve_Succeeded || status(model)==:Optimal || status(model)==:NearOptimal || status(model)==:KnitroError || ( status(model)==:UserLimit && solver !=:Mosek ) || status(model)==:FeasibleApproximate 
        q = Vector{Float64}( getsolution(model) )
        grad = zeros(Float64, myeval.n)
        InfDecomp.∇f(myeval, grad, q, Float64(.0))

        lambda = Vector{Float64}( getconstrduals(model) )

        mu = grad - myeval.Gt*lambda

        fstat.mu_nonneg_viol = -minimum(mu)

        fstat.complementarity_max = maximum( abs.(mu) .* abs.(q) )
        fstat.complementarity_sum = sum( abs.(mu) .* abs.(q) )
    else
        q = Vector{Float64}( ill_sol )
        fstat.mu_nonneg_viol = -10
        fstat.complementarity_max = -10
        fstat.complementarity_sum = -10
    end# ^if status

    fstat.obj_val = eval_f(myeval,q,Float64(0))

    fstat.q_nonnegativity = -minimum(q)

    equation = (q'*myeval.Gt)' - myeval.rhs
    fstat.marginals_1   = norm(equation,1)
    fstat.marginals_2   = norm(equation,2)
    fstat.marginals_Inf = norm(equation,Inf)
    # fstat.entropy_X   = entropy_X(myeval,q,Float64(0))
    # fstat.MI_X_YZ     = fstat.entropy_X + fstat.obj_val
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

    fstat.CI = -(eval_f(myeval,p,Float64(0)) + fstat.obj_val)/log(2)
    fstat.SI = SI(myeval,q,Float64(0))
    fstat.UI_Y = I_X_Y__Z(myeval,q,Float64(0))
    fstat.UI_Z = I_X_Z__Y(myeval,q,Float64(0))
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
    ;
end #^ check_feasibility()

end #^ module

; # EOF

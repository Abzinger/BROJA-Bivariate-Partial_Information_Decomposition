# infdecomp.jl
module InfDecomp
export decomp

using MathProgBase
# http://mathprogbasejl.readthedocs.io/en/latest/nlp.html

type My_Eval <: MathProgBase.AbstractNLPEvaluator
    n_x     :: Int32
    n_y     :: Int32
    n_z     :: Int32

    n     :: Int32         # number of variables
    m     :: Int32         # number of marginal equations

    varidx  :: Array{Int32,3}                   # 0 if variable not present; otherwise idx of var in {1,...,n}
    xyz     :: Vector{Tuple{Int32,Int32,Int32}} # xyz-triple of varidx

    eqidx   :: Dict{ Tuple{String,Int32,Int32},  Int32} # first idx is "xy" or "xz". Index of eqn in {1,...,m}
    mr_eq   :: Vector{ Tuple{String,Int32,Int32} }      # ("xy",x,y) / ("yz",y,z) of an eqn

    marg_xy :: Array{Float64,2}
    marg_xz :: Array{Float64,2}
    marg_x  :: Array{Float64,1}

    rhs   :: Vector{Float64}                # m-vector
    Gt    :: SparseMatrixCSC{Float64,Int32} # G^T n x m; transpose of constraint matrix
    Gt_K  :: Vector{Int32}                  # findn()-info of G^T: rows
    Gt_L  :: Vector{Int32}                  # findn()-info of G^T: columns

    TmpFloat       :: DataType
    bigfloat_nbits :: Int32
end

function create_My_Eval(q::Array{Float64,3})
    const n_x::Int32 = size(q,1);
    const n_y::Int32 = size(q,2);
    const n_z::Int32 = size(q,3);

    # Create marginals
    marg_xy::Array{Float64,2} = zeros(n_x,n_y)
    marg_xz::Array{Float64,2} = zeros(n_x,n_z)
    marg_x::Array{Float64,1}  = zeros(n_x)
    for x in 1:n_x
        for y in 1:n_y
            for z in 1:n_z
                marg_xy[x,y] += q[x,y,z]
                marg_xz[x,z] += q[x,y,z]
                marg_x[x]    += q[x,y,z]
            end
        end
    end

    # Find the variables
    varidx::Array{Int32,3}                   = zeros(Bool,size(q));
    xyz   ::Vector{Tuple{Int32,Int32,Int32}} = [ (0,0,0) for i in 1:n_x*n_y*n_z ]
    n::Int32 = 0
    for x in 1:n_x
        for y in 1:n_y
            for z in 1:n_z
                if marg_xy[x,y] > 0  &&  marg_xz[x,z] > 0
                    n += 1
                    varidx[x,y,z] = n
                    xyz[n]        = (x,y,z)
                else
                    varidx[x,y,z] = 0
                end#if
            end
        end
    end

    # Find the equations
    eqidx = Dict{ Tuple{String,Int32,Int32},Int32}() # first idx is "xy" or "xz"
    mr_eq ::Vector{ Tuple{String,Int32,Int32} }   = [ ("",0,0)   for i in 1:n_x*(n_y+n_z) ]
    m::Int32 = 0
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

    Gt::SparseMatrixCSC{Float64,Int32} = sparse(denseGt)
    local Gt_K::Array{Int32,1}
    local Gt_L::Array{Int32,1}
    (Gt_K,Gt_L) = findn(Gt)

    TmpFloat       :: DataType  = BigFloat
    bigfloat_nbits :: Int32     = 256


    return My_Eval(n_x,n_y,n_z, n,m, varidx,xyz, eqidx,mr_eq, marg_xy,marg_xz,marg_x, rhs, Gt,Gt_K,Gt_L,  TmpFloat,bigfloat_nbits)
    ;
end


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
# import MathProgBase.getreducedcosts
# import MathProgBase.getconstrduals

# ------------
# B a s i c s
# ------------

# initialize()
function initialize(e::My_Eval, requested_features::Vector{Symbol})
    for feat in requested_features
        if feat ∉ e.features_list
            error("infdecomp.jl: initialize():\n-   JuliaOpt:MathProgBase is asking for a feature ($feat) that I don't have.-   Maybe use another solver?")
        end
    end
end


# features_available()
features_available(::My_Eval) = features_list

# Properties:
isobjlinear(::My_Eval)           = false
isobjquadratic(::My_Eval)        = false
isconstrlinear(::My_Eval, ::Any) = true


# ------------------------------------------
# E v a l u a t i o n :   0 t h   o r d e r
# ------------------------------------------
function condEntropy{TFloat}(e::My_Eval, p::Vector{Float64}, dummy::TFloat)   :: TFloat
    # m_yz = marg_yz(e,p,dummy)
    s::TFloat = TFloat(0.)
    for y = 1:e.n_y
        for z = 1:e.n_z
            P_yz = TFloat(0.)
            for x = 1:e.n_x
                i = e.varidx(x,y,z)
                if i>0
                    P_yz += e.p[i]
                end
            end
            for x = 1:e.n_x
                i = e.varidx(x,y,z)
                if i>0
                    p_xyz = e.p[i]
                    s  +=  (  p_xyz ≤ 0 ?   TFloat(0.)   :   - p_xyz*log( p_xyz / P_yz )   )
                end
            end
        end
    end
    return s
    ;
end

function eval_f(e::My_Eval, x::Vector{Float64}) :: Float64
    local condent::Float64
    if TmpFloat==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        condent = condEntropy(e,x,BigFloat())
        setprecision(BigFloat,prc)
    else
        condent = condEntropy(e,x,Float64())
    end

    return -condent
    ;
end

# eval_g --- eval of constraint into g
function eval_g(e::My_Eval, g::Vector{Float64}, x::Vector{Float64})  :: Void
    g .= reshape( reshape(x,1,e.n)*e.Gt , e.m, ) .- rhs
    return nothing
    ;
end # eval_g()


# ------------------------------------------
# E v a l u a t i o n :   1 s t   o r d e r
# ------------------------------------------

function ∇f{TFloat}(e::My_Eval, grad::Vector{Float64}, p::Vector{Float64}, dummy::TFloat) :: Void
    for y = 1:e.n_y
        for z = 1:e.n_z
            # make marginal P(*yz)
            P_yz::TFloat = TFloat(0.)
            for x = 1:e.n_x
                i = e.varidx(x,y,z)
                if i>0
                    P_yz += e.p[i]
                end
            end
            # make log-expressions  log( P(xyz) / P(*yz) )
            for x = 1:e.n_x
                i = e.varidx(x,y,z)
                if i>0
                    P_xyz::TFloat = TFloat( e.p[i] )
                    grad[i] = Float64(  P_xyz ≤ 0 ?  TFloat(0.)  : log( P_xyz / P_yz )  )
                end
            end
        end# for y
    end# for x
    ;
end # ∇f()

# eval_grad_f --- eval gradient of objective function
function eval_grad_f(e::My_Eval, g::Vector{Float64}, x::Vector{Float64}) :: Void
    if TmpFloat==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        ∇f(e,g,x,BigFloat())
        setprecision(BigFloat,prc)
    else
        ∇f(e,g,x,Float64())
    end
    return nothing
    ;
end # eval_grad_f()


# Constraint Jacobian
# jac_structure() --- zero-nonzero pattern of constraint Jacobian
jac_structure(e::My_Eval) :: Tuple(Vector{Int32},Vector{Int32})   =  ( e.Gt_L , e.Gt_K )
# Note: Gt is transposed, so K and L are swapped [K: rows of G^T; L: columns of G^T]

# eval_jac_g() --- constraint Jacobian   -> J
function eval_jac_g(e::My_Eval, J::Vector{Float64}, x::Vector{Float64}) :: Void
    J .= e.Gt.nzval
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
function hesslag_structure(e::My_Eval)  :: Tuple(Vector{Int32},Vector{Int32})
    counter = 0
    for y = 1:e.n_y
        for z = 1:e.n_z
            # Start with the diagonal
            for x = 1:e.n_x
                i = e.varidx[x,y,z]
                if i>0
                    counter += 1
                    K[counter] = i
                    L[counter] = i
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
                        K[counter] = i_x
                        L[counter] = i_u
                    end
                end
            end#for z
        end# for y

    end #^ for yz
    return nothing
    ;
end # hesslag_structure()


function Hess{TFloat}(e::My_Eval, H::Vector{Float64}, x::Vector{Float64}, σ::Float64, dummy::TFloat) :: Void
    counter = 0
    for y = 1:e.n_y
        for z = 1:e.n_z
            # make marginal P(*yz)
            P_yz  ::TFloat = TFloat(0.)
            for x = 1:n_X
                i = e.varidx[x,y,z]
                if i>0
                    P_yz += e.p[i]
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
                    P_xyz = e.p[i]
                    H[counter] = Float64(  TFloat(σ)*( P_yz - P_xyz )/(  P_yz * P_xyz )  )
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
    if TmpFloat==BigFloat
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        Hess(e,H,x,σ,BigFloat())
        setprecision(BigFloat,prc)
    else
        Hess(e,H,x,σ,Float64())
    end
    return nothing
    ;
end # eval_hesslag()


# eval_hesslag() --- ( Hessian [wrt x] of the Lagrangian ) * v
# function eval_hesslag_prod{T<:AbstractFloat}(e::My_Eval, h::Vector{T}, x::Vector{T}, v::Vector{T}, σ::T, μ::Vector{T}) :: Void
#     h .= σ .* Hf(x)*v
#     return nothing
# end

# ----------------
# U S I N G   I T
# ----------------

numo_vars(e::My_Eval)                   :: Int32            = e.n
numo_constraints(e::My_Eval)            :: Int32            = e.m
constraints_lowerbounds_vec(e::My_Eval) :: Vector{Float64}  = zeros{Float64}(e.m)
constraints_upperbounds_vec(e::My_Eval) :: Vector{Float64}  = zeros{Float64}(e.m)
vars_lowerbounds_vec(e::My_Eval)        :: Vector{Float64}  = zeros{Float64}(e.n)
vars_upperbounds_vec(e::My_Eval)        :: Vector{Float64}  = [Inf  for j=1:e.n]


type Set_Data{T1,T2,T3}
    Xidx :: Dict{T1,Int32}
    Yidx :: Dict{T2,Int32}
    Zidx :: Dict{T3,Int32}
    X    :: Vector{T1}
    Y    :: Vector{T1}
    Z    :: Vector{T1}
end

function create_Set_Data{T1,T2,T3}(pdf::Dict{Tuple{T1,T2,T3},Float64}) :: Set_Data{T1,T2,T3}
    Xidx = Dict{T1,Int32}()
    Yidx = Dict{T2,Int32}()
    Zidx = Dict{T3,Int32}()

    for (xyz,val) in pdf
        x,y,z = xyz
        if x ∉ keys(Xidx)
            Xidx[x] = length(Xidx)+1
        end
        if y ∉ keys(Yidx)
            Yidx[y] = length(Yidx)+1
        end
        if z ∉ keys(Zidx)
            Zidx[z] = length(Zidx)+1
        end
    end

    X ::Vector{T1} = [ Xidx[1]  for i = 1:length(Xidx) ]
    Y ::Vector{T1} = [ Yidx[1]  for i = 1:length(Yidx) ]
    Z ::Vector{T1} = [ Zidx[1]  for i = 1:length(Zidx) ]

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
end


function create_stuff{T1,T2,T3}(pdf::Dict{Tuple{T1,T2,T3},Float64}, tmpFloat::DataType=BigFloat) :: Tuple{ Set_Data{T1,T2,T3}, My_Eval }
    sd = create_Set_Data(pdf)

    Q = zeros(Float64,length(X),length(Y),length(Z))
    for (xyz,val) in pdf
        Q[ sd.X[xyz[1]], sd.Y[xyz[2]], sd.Z[xyz[3]] ] = val
    end

    e = create_My_Eval(Q)

    return (sd,e)
    ;
end #^ create_stuff

using Base.Test


function do_it{T1,T2,T3}(pdf::Dict{Tuple{T1,T2,T3},Float64}, tmpFloat::DataType=BigFloat, solver=MathProgBase.defaultNLPsolver)
    const model = MathProgBase.NonlinearModel(solver)
    const sd,myeval = create_stuff(pdf)
    const lb = myeval.constraints_lowerbounds_vec
    const ub = myeval.constraints_upperbounds_vec
    const l = myeval.vars_lowerbounds_vec
    const u = myeval.vars_upperbounds_vec

    MathProgBase.loadproblem!(model, myeval.n, myeval.m, l, u, lb, ub, :Max, myeval)

    MathProgBase.optimize!(model)
    stat = MathProgBase.status(model)


    # @test stat == :Optimal
    # @show MathProgBase.getobjval(model)
    # objvaldist = abs( log(n) - MathProgBase.getobjval(model) )*1.e-9
    # println("Distance from true optimum (in 1.e-9): $objvaldist")

    # x = MathProgBase.getsolution(model)
    # dist::Float64 = 0.
    # for j=1:n
    #     dist += abs2(x[j]-1./n)
    # end
    # dist = sqrt(dist)*1.e-9
    # println("Norm Distance from true optimal value (in 1.e-9): $dist")

    # # @test_approx_eq_eps MathProgBase.getobjval(model) log(n)  1.e-300

#    # Test that a second call to optimize! works
#    MathProgBase.setwarmstart!(m,[1,5,5,1])
#    MathProgBase.optimize!(m)
#    stat = MathProgBase.status(m)
#    @test stat == :Optimal
end

end # module

; # EOF

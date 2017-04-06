# infdecomp.jl
module InfDecomp
export decomp, verify

using MathProgBase

type My_Eval <: MathProgBase.AbstractNLPEvaluator
    n_x     :: Int32
    n_y     :: Int32
    n_z     :: Int32

    varidx  :: Array{Int32,3} # 0 if variable not present; otherwise idx of var
    xyz     :: Vector{Tuple{Int32,Int32,Int32}} # xyz-triple of varidx

    eqidx   :: Dict{ Tuple{String,Int32,Int32},Int32} # first idx is "xy" or "xz"
    mr_eq   :: Vector{ Tuple{String,Int32,Int32} }    # ("xy",x,y) / ("yz",y,z) of an eqn

    marg_xy :: Array{Float64,2}
    marg_xz :: Array{Float64,2}
    marg_x  :: Array{Float64,1}

    n     :: Int32         # number of variables
    m     :: Int32         # number of marginal equations


    rhs   :: Vector{Float64} # m-vector
    Gt    :: SparseMatrixCSC{Float64,Int32} # G^T n x m; transpose of constraint matrix
    Gt_K  :: Vector{Int32} # findn()-info of G^T: rows
    Gt_L  :: Vector{Int32} # findn()-info of G^T: columns

    TmpFloat       :: DataType
    bigfloat_nbits :: Int32
end

function create_My_Eval(q::Array{Float64,3})
    if ndims(q) != 3
        print("Need 3 dimensions in q\n");
        return;
    end
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
    varidx::Array{Int32,3} = zeros(Bool,size(q));
    xyz   :: Vector{Tuple{Int32,Int32,Int32}} = [ (0,0,0) for i in 1:n_x*n_y*n_z ]
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
    mr_eq::Vector{ Tuple{String,Int32,Int32} }    = [ ("",0,0)   for i in 1:n_x*(n_y+n_z) ]
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


    return My_Eval(n_x,n_y,n_z, varidx,xyz, eqidx,mr_eq, marg_xy,marg_xz,marg_x, n,m, rhs, Gt,Gt_K,Gt_L,  TmpFloat,bigfloat_nbits)
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

function condEntropy{TFloat}(e::My_Eval, x::Vector{Float64}, dummy::TFloat)   :: TFloat
    P = reshape(x,sz_X,sz_YZ)
    s::TFloat  = TFloat(0.)
    for yz = 1:sz_YZ
        # make marginal P(*yz)
        P_yz::TFloat = TFloat(0.)
        for x = 1:sz_X
            P_yz += P[x,yz]
        end
        # make and add log-expressions  P(xyz)* log( P(xyz) / P(*yz) )
        for x = 1:sz_X
            P_xyz::TFloat = TFloat( P[x,yz] )
            P_xyz ≤ 0 || (   s += P_xyz * log( P_xyz / P_yz )   )
        end
    end
    return s
    ;
end

function eval_f(e::My_Eval, x::Vector{Float64}) :: Float64
    local retval::Float64
    if (TmpFloat==BigFloat)
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        retval = condEntropy(e,x,BigFloat())
        setprecision(BigFloat,prc)
    else
        retval = condEntropy(e,x,Float64())
    end

    return retval
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

function ∇f{TFloat}(e::My_Eval, g::Vector{Float64}, x::Vector{Float64}, dummy::TFloat) :: Void
    grad = reshape(g, sz_X,sz_YZ)
    P    = reshape(x, sz_X,sz_YZ)
    for yz = 1:sz_YZ
        # make marginal P(*yz)
        P_yz::TFloat = TFloat(0.)
        for x = 1:sz_X
            P_yz += P[x,yz]
        end
        # make log-expressions  log( P(xyz) / P(*yz) )
        for x = 1:sz_X
            P_xyz::TFloat = TFloat( P[x,yz] )
            grad[x,yz] = Float64(  P_xyz ≤ 0 ?  TFloat(0.)  : log( P_xyz / P_yz )  )
        end
    end
    ;
end # ∇f()

# eval_grad_f --- eval gradient of objective function
function eval_grad_f(e::My_Eval, g::Vector{Float64}, x::Vector{Float64}) :: Void
    if (TmpFloat==BigFloat)
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
jac_structure(e::My_Eval) ::Tuple(Vector{Int32},Vector{Int32})  = ( e.Gt_L , e.Gt_K ) # Note: Gt is transposed, so K and L are swapped


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
    P   = reshape(x, sz_X,sz_YZ)
    counter = 1
    for yz = 1:sz_YZ
        # Start with the diagonal
        for x = 1:sz_X
            K[counter] = varidx(e,x,yz)
            L[counter] = varidx(e,x,yz)
            counter += 1
        end
        # Now off-diagonal.
        # H is  treated as a symmetric matrix, though (*): if both (k,l) and (l,k) are present, their values will be added!
        for x = 1:sz_X
            for u = 1:sz_X
                if x ≠ u
                    K[counter] = varidx(e,x,yz)
                    L[counter] = varidx(e,u,yz)
                    counter += 1
                end
            end
        end

    end #^ for yz
    return nothing
    ;
end # hesslag_structure()


function Hess{TFloat}(e::My_Eval, H::Vector{Float64}, x::Vector{Float64}, σ::Float64, dummy::TFloat) :: Void
    P    = reshape(x, sz_X,sz_YZ)
    counter = 1
    for yz = 1:sz_YZ
        # make marginal P(*yz)
        P_yz::TFloat = TFloat(0.)
        for x = 1:sz_X
            P_yz += P[x,yz]
        end
        # now: for all pairs x,u we have:
        # if x ≠  u:   -1/P_yz;
        # if x == u:   ( P_yz - P(xyz) )/(  P_yz * P(xyz) )

        # Start with the diagonal
        for x = 1:sz_x
            P_xyz::TFloat = TFloat( P[x,yz] )
            if P_xyz > 0
                H[counter] = Float64(  σ*( P_yz - P_xyz )/(  P_yz * P_xyz )  )
            else
                H[counter] = 0. # always the same question: does this make sense?!?
            end
            counter += 1
        end
        # Now off-diagonal.
        # H must be stored as a symmetric matrix, though, so at most one of (k,l) , (l,k) is present.
        for x = 1:sz_X
            for u = 1:sz_X
                if x ≠ u
                    H[counter] =  ½ * Float64(  -σ/P_yz  ) # the ½ is because the matrix is treated as symmetric, see (*) above.
                    counter += 1
                end
            end
        end

    end #^ for yz
    return nothing
    ;
end # eval_hesslag()



# eval_hesslag() --- Hessian [wrt x] of the Lagrangian
function eval_hesslag(e::My_Eval, H::Vector{Float64}, x::Vector{Float64}, σ::Float64, μ::Vector{Float64}) :: Void
    if (TmpFloat==BigFloat)
        prc = precision(BigFloat)
        setprecision(BigFloat,e.bigfloat_nbits)
        Hess(e,H,x,σ,BigFloat())
        setprecision(BigFloat,prc)
    else
        Hess(e,H,x,σ,BigFloat())
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


function create_stuff{T1,T2,T3,FLOAT<:AbstractFloat}(pdf::Dict{Tuple{T1,T2,T3},FLOAT}, tmpFloat::DataType=Float64) :: Tuple{Set_Data{T1,T2,T3},My_Eval{FLOAT}}
    e = My_Eval{FLOAT}()
    d = Set_Data{T1,T2,T3}()

    for xyz ∈ keys(pdf)
        ( pdf[xyz] > 0 )||(  delete!(pdf,xyz)  )
    end


    @assert ( length(d.X) ≥ 2 ) "|Range(X)| ≥ 2 needed"
    @assert ( length(d.X) ≥ 2 ) "|Range(Y)| ≥ 2 needed"
    @assert ( length(d.X) ≥ 2 ) "|Range(Z)| ≥ 2 needed"

    e.n     = length(pdf)
    e.m     = length(d.X)*length(d.Y) + length(d.X)*length(d.Z)
    e.sz_X  = length(d.X)
    e.sz_YZ = length(d.Y)*length(d.Z)

    e.rhs = zeros(FLOAT,e.m)

    # Make tmpG and rhs
    tmpG = zeros(Float64,e.m,e.n)
    begin
        local row::Int = 1
        begin # xy-marginals
            for (i_x,x) in enumerate(d.X)
                for (i_y,y) in enumerate(d.Y)
                    for (i_z,z) in enumerate(d.Z)
                        tmpG[row, varidx(e,i_x,length(d.Z)*(i_y-1)+i_z)] = true
                        e.rhs[row] += get(pdf,(x,y,z),0.)
                    end
                    row += 1
                end
            end
        end #^ xy-marginals
        begin # xz-marginals
            for (i_x,x) in enumerate(d.X)
                for (i_z,z) in enumerate(d.Z)
                    for (i_y,y) in enumerate(d.Y)
                        tmpG[row, varidx(e,i_x,length(d.Z)*(i_y-1)+i_z)] = true
                        e.rhs[row] += get(pdf,(x,y,z),0.)
                    end
                    row += 1
                end
            end
        end #^ xz-marginals
    end #^ begin

    e.Gt = sparse(tmpG')
    e.Gt_K,e.Gt_L = findn(e.Gt)

    tmpG = nothing

    return (d,e)
    ;
end #^ create_stuff

using Base.Test

function mytest(n::Int32, solver=MathProgBase.defaultNLPsolver)

    const model = MathProgBase.NonlinearModel(solver)
    const myeval = create_My_Eval(n)
    const lb = myeval.constraints_lowerbounds_vec
    const ub = myeval.constraints_upperbounds_vec
    const l = myeval.vars_lowerbounds_vec
    const u = myeval.vars_upperbounds_vec

    MathProgBase.loadproblem!(model, myeval.numo_vars, myeval.numo_constraints, l, u, lb, ub, :Max, myeval)

    MathProgBase.optimize!(model)
    stat = MathProgBase.status(model)

    @test stat == :Optimal
    @show MathProgBase.getobjval(model)
    objvaldist = abs( log(n) - MathProgBase.getobjval(model) )*1.e-9
    println("Distance from true optimum (in 1.e-9): $objvaldist")

    x = MathProgBase.getsolution(model)
    dist::Float64 = 0.
    for j=1:n
        dist += abs2(x[j]-1./n)
    end
    dist = sqrt(dist)*1.e-9
    println("Norm Distance from true optimal value (in 1.e-9): $dist")

    # @test_approx_eq_eps MathProgBase.getobjval(model) log(n)  1.e-300

#    # Test that a second call to optimize! works
#    MathProgBase.setwarmstart!(m,[1,5,5,1])
#    MathProgBase.optimize!(m)
#    stat = MathProgBase.status(m)
#    @test stat == :Optimal
end

end # module

; # EOF

module InfDecomp_Base

export My_Eval, condEntropy, create_My_Eval, Set_Data, create_stuff

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
    Gt_K  :: Vector{Int64}                  # findn()-info of G^T: rows (vars)
    Gt_L  :: Vector{Int64}                  # findn()-info of G^T: columns (equations)

    TmpFloat       :: DataType
    bigfloat_nbits :: Int64
end

function condEntropy{TFloat_2,TFloat}(e::My_Eval, p::Vector{TFloat_2}, dummy::TFloat)   :: TFloat
    # m_yz = marg_yz(e,p,dummy)
    s::TFloat = TFloat(0.)
    for y = 1:e.n_y
        for z = 1:e.n_z
            P_yz = TFloat(0.)
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
                    s  +=  (  (p_xyz ≤ 0 || P_yz  ≤ 0)  ?   TFloat(0.)   :   - p_xyz*log( p_xyz / P_yz )   )
                end
            end
        end
    end
    return s
    ;
end


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


end #^ module InfDecomp_base
;#EOF

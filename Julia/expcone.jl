# expcone.jl -- part of the Information Decomposition project
module ExpCone

#incude("infdecomp_base.jl")
using InfDecomp_Base
using MathProgBase
#using CPUTime

type Solution_Stats
    optimum     :: Float64
    primal_val  :: Float64
#    primal_feas :: Float64
#    dual_val    :: Float64
#    dual_feas   :: Float64
    q           :: Vector{Float64}
    time        :: Float64
    model       :: MathProgBase.AbstractMathProgModel
end

function model_L(D::InfDecomp_Base.My_Eval, solver::MathProgBase.AbstractMathProgSolver) :: Solution_Stats
    # min c^T v     #  min (-1,0,0)^T (r^T,q^T,p^T)
    # s.t.          #  s.t.
    # b - Av = 0    #  q_{*xy} = p_{xyz} ;  + marginal eqns
    # v ∈ EXP       #  (r,q,p) ∈ EXP

    local n_vars::Int = D.n*3            # r_{xyz}, q_{xyz}, p_{xyz} ;  ∀ xyz ∈ Variables
                                         # r=3*i-2 ; q=3*i-1 ; p=3*i
    local n_cons::Int = D.n + D.m        #    0          q_{*yz} - p_{xyz}  = 0     ; ∀ xyz
                                         # +  b^y_{xy} - q_{xy*}            = 0     ; ∀ xy
                                         # +  b^z_{xz} - q_{x*z}            = 0     ; ∀ xz

    function make_c_A_b() :: Tuple{Vector{Float64},SparseMatrixCSC{Float64},Vector{Float64}}
        local c::Vector{Float64} = zeros(n_vars)
        local b::Vector{Float64} = zeros(n_cons)

        # fill c: -1 on the r_{xyz} vars:
        for xyz = 1:D.n
            c[3*xyz-2] = -1.
        end

        # fill b:  (0, b^y, b^z)
        for i = 1:D.m
            b[ D.n + i ]  = D.rhs[i]
        end

        #
        # fill A
        # ------

        Eqn   = Array{Int,1}()
        Var   = Array{Int,1}()
        Coeff = Array{Float64,1}()

        # q=p equations:
        for xyz = 1:D.n
            eqn = xyz

            p_xyz = 3*xyz
            append!(Eqn,   eqn)
            append!(Var,   p_xyz)
            append!(Coeff, -1.)

            (x,y,z) = D.xyz[xyz]
            for uvw = 1:D.n
                (u,v,w) = D.xyz[uvw]
                if (v,w)==(y,z)
                    q_uvw = 3*uvw-1
                    append!(Eqn,   eqn)
                    append!(Var,   q_uvw)
                    append!(Coeff, 1.)
                end
            end
        end

        # q_{xy*} = b^y_{xy}
        (K,L,V) = findnz(D.Gt)
        for eq = 1:length(K)
            append!(Eqn,   D.n+L[eq])
            append!(Var,   3*K[eq]-1)
            append!(Coeff, V[eq])
        end
        A = sparse(Eqn,Var,Coeff)
        @show A
        return c,A,b
    end #^ make_c_b_A()


    c,A,b = make_c_A_b()

    constr_cones = [ ( :Zero, 1:n_cons ) ]
    var_cones = []
    for xyz = 1:D.n
        rqp_xyz = [3*xyz-2,3*xyz-1,3*xyz]
        append!(var_cones, [(:ExpPrimal, rqp_xyz)] )
    end


    m = MathProgBase.ConicModel(solver)
    MathProgBase.loadproblem!(m, c,A,b, constr_cones, var_cones)

    # CPUtic()                        # TIC

    MathProgBase.optimize!(m)

    # stats_time = CPUtoq()           # TOQ


    @show MathProgBase.status(m)
    @show stats_time = MathProgBase.getsolvetime(m)


    if MathProgBase.status(m) == :Optimal || MathProgBase.status(m) == :Suboptimal || MathProgBase.status(m) == :Error
        q = Vector{Float64}(D.n)
        for xyz in 1:D.n
            q[xyz] = MathProgBase.getsolution(m)[ 3*xyz-1 ]
        end

        eval_val::Float64 = InfDecomp_Base.condEntropy(D,q,big(0.))

        @show eval_val, MathProgBase.getobjval(m)
        if abs(eval_val + MathProgBase.getobjval(m)) > 1.e-3
            println("BIG DIFFERENCE IN SOLUTION VALUE!! :-( :-( :-( :-( :-( :-(\nProbably a BUG!!!!!!!!!!!!!!!!!!!!!")
        end

        stats = Solution_Stats(-eval_val,MathProgBase.getobjval(m),q,stats_time,m)
        return stats
    else
        print("Fuck!!")
        stats = Solution_Stats(999,999,[],stats_time,m)
        return stats
    end
    ; #
end #^ model_L


######################################################################################################################################################
#
######################################################################################################################################################


function model_S(D::InfDecomp_Base.My_Eval, solver::MathProgBase.AbstractMathProgSolver) :: Solution_Stats
    #    min (-1,0,0)(r,q,p)
    # s.t.
    #    q_{*yz} - p_{yz}         = 0     p/q coupling
    #    b - Gq                   = 0     marginal eqns
    #    (r_{xyz},q_{xyz},p_{yz}) ∈ EXP   linear EXPCONE_ieq
    #
    #    rqp free variables

    # count yz's
    local idx_of_yz = Dict{Tuple{Int64,Int64}, Int64}()

    local n_yz = 0
    for xyz = 1:D.n
        (x,y,z) = D.xyz[xyz]
        if (y,z) ∉ keys(idx_of_yz)
            n_yz += 1
            idx_of_yz[(y,z)] = n_yz
        end
    end

    local yz_of_idx = Vector{Tuple{Int64,Int64}}(n_yz)

    for (yz,i) in idx_of_yz
        yz_of_idx[i] = yz
    end


    local n_vars::Int = D.n + D.n + n_yz      # r_{xyz} ∀xyz;  q_{xyz} ∀xyz;  p_{yz} ∀yz
                                              # r=1...n; q=n+1...2n; p=2n+1...2n+n_yz
    local n_cons::Int = n_yz + D.m + 3*D.n    #   0           - ( q_{*yz}  -p_{yz}        ) ∈ {0}   ; ∀ yz
                                              #   / b^y_{xy}  - ( q_{xy*}                 ) ∈ {0}   ; ∀ xy  \  | D.m
                                              #   \ b^z_{xz}  - ( q_{x*z}                 ) ∈ {0}   ; ∀ xz  /  |
                                              #   0           - (-r_{xyz},-q_{xyz},-p_{yz}) ∈ EXPC  ; ∀ xyz

    function make_c_A_b() :: Tuple{Vector{Float64},SparseMatrixCSC{Float64},Vector{Float64}}
        local c::Vector{Float64} = zeros(n_vars)
        local b::Vector{Float64} = zeros(n_cons)

        # fill c: -1 on the r_{xyz} vars:
        for xyz = 1:D.n
            c[xyz] = -1.
        end

        # fill b:  (0, b^y b^z, 0)
        for i = 1:D.m
            b[ n_yz + i ]  = D.rhs[i]
        end

        #
        # fill A
        # ------

        Constr = Array{Int,1}()
        Var    = Array{Int,1}()
        Coeff  = Array{Float64,1}()

        # q=p equations:
        for yz = 1:n_yz
            eqn = yz
            p_yz = 2*D.n + yz

            append!(Constr, eqn)
            append!(Var,    p_yz)
            append!(Coeff,  -1.)

            (y,z) = yz_of_idx[yz]

            for uvw = 1:D.n
                (u,v,w) = D.xyz[uvw]
                if (v,w) == (y,z)
                    q_uvw = D.n + uvw
                    append!(Constr, eqn)
                    append!(Var,    q_uvw)
                    append!(Coeff,  1.)
                end
            end
        end

        # q_{xy*} = b^y_{xy}
        (K,L,V) = findnz(D.Gt)
        for eq = 1:D.m
            append!(Constr, n_yz + L[eq])
            append!(Var,    D.n + K[eq])
            append!(Coeff,  V[eq])
        end

        # EXPC ieqs:
        for xyz = 1:D.n
            cns_base = n_yz + D.m + 3*xyz-3
            r_xyz = xyz
            q_xyz = D.n + xyz
            (x,y,z) = D.xyz[xyz]
            p_yz = 2*D.n + idx_of_yz[ (y,z) ]

            append!(Constr, cns_base +1)
            append!(Var,    r_xyz)
            append!(Coeff,  -1.)

            append!(Constr, cns_base +2)
            append!(Var,    q_xyz)
            append!(Coeff,  -1.)

            append!(Constr, cns_base +3)
            append!(Var,    p_yz)
            append!(Coeff,  -1.)
        end

        A = sparse(Constr,Var,Coeff)

        return c,A,b
    end #^ make_c_b_A()

    (c,A,b) = make_c_A_b()

    var_cones = [ ( :Free, 1:n_vars) ]

    constr_cones = []
    append!(constr_cones, [ ( :Zero, 1:n_yz             ) ] )
    append!(constr_cones, [ ( :NonNeg, (n_yz+1):(n_yz+D.m)) ] )
    for xyz = 1:D.n
        rqp = [ n_yz+D.m+3*xyz-2, n_yz+D.m+3*xyz-1, n_yz+D.m+3*xyz ]
        append!(constr_cones, [(:ExpPrimal, rqp)] )
    end

    m = MathProgBase.ConicModel(solver)
    MathProgBase.loadproblem!(m, c,A,b, constr_cones, var_cones)


    MathProgBase.optimize!(m)


    @show MathProgBase.status(m)

    if MathProgBase.status(m) == :Optimal
        q = Vector{Float64}(D.n)
        for xyz in 1:D.n
            q[xyz] = MathProgBase.getsolution(m)[ D.n + xyz ]
        end
        eval_val::Float64 = InfDecomp_Base.condEntropy(D,q,big(0.))

        @show eval_val, MathProgBase.getobjval(m)
        if abs(eval_val + MathProgBase.getobjval(m)) > 1.e-3
            println("BIG DIFFERENCE IN SOLUTION VALUE!! :-( :-( :-( :-( :-( :-(\nProbably a BUG!!!!!!!!!!!!!!!!!!!!!")
        end

        stats = Solution_Stats(-eval_val,MathProgBase.getobjval(m),q,m)

        return stats
    else
        print("Fuck!!")
        stats = Solution_Stats(0,0,[],m)
        return stats
    end
    ; #
end #^ model_S

################################################################################################################################################
#
###############################################################################################################################################

function model_D(D::InfDecomp_Base.My_Eval, solver::MathProgBase.AbstractMathProgSolver) :: Solution_Stats
    # min c^T v     #  min (0,b)^T (u^T,l^T)
    # s.t.          #  s.t.
    #               #     
    # f(v) ∈ EXP*   #  (-1,f(u,l),-u) ∈ EXP*

    local n_vars::Int = D.n + D.m        # u_{xyz} ∀ xyz; l_{xy} ∀ xy; l_{xz}  ∀ xz;
                                         # u=i                  ;      l=D.n+i
    local n_cons::Int = D.n*3            # u_{*yz} + ln(-u_{xyz}) + l_{xy} + l_{xz} \le -1     ; ∀ xyz

    function make_c_A_b() :: Tuple{Vector{Float64},SparseMatrixCSC{Float64},Vector{Float64}}
        local c::Vector{Float64} = zeros(n_vars)
        local b::Vector{Float64} = zeros(n_cons)

        println("the elements are ", D.xyz)
        # fill c: (0, b^y, b^z):
        for i = 1:D.m
            c[ D.n + i ]  = D.rhs[i]
        end

        # fill b: (-1,0,0):
        for xyz = 1:D.n
            b[xyz] = -1.
        end

        #
        # fill A
        # ------

        Constr   = Array{Int,1}()
        Var   = Array{Int,1}()
        Coeff = Array{Float64,1}()

        # -1  EXP* ieqs:
        for xyz = 1:D.n
            ieqn = xyz
            u_xyz = xyz
            append!(Constr, ieqn)
            append!(Var, u_xyz)
            append!(Coeff, 0.)
        end

        # u_{*yz} + l_xy + l_xz EXP* ieqs:
        for xyz = 1:D.n
            (x,y,z) = D.xyz[xyz]
            ieqn = xyz
            # add u_{*yz}
            for rpq = 1:D.n
                (r,p,q) = D.xyz[rpq]
                if (p,q) == (y,z)
                    u_rpq = rpq
                    println("u_rpq = ($r,$p,$q) with idx= $u_rpq")
                    append!(Constr, D.n + ieqn)
                    append!(Var,    u_rpq)
                    append!(Coeff,  -1.)
                end
            end
            # add l_xy
            if D.marg_xy[x,y] > 0
                append!(Constr, D.n + ieqn)
                println("l_xy = ($x ,$y) with idx= ", D.n + D.eqidx["xy",x,y])
                append!(Var,    D.n + D.eqidx["xy",x,y])
                append!(Coeff,  -1.)
            end
            # add l_xz
            # !! D.eqidx adds xz after xy. So no need
            # to have D.n + mr_xy + D.eqidx["xz",x,z] !!
            if D.marg_xz[x,z] > 0
                append!(Constr, D.n + ieqn)
                println("l_xz = ($x ,$z) with idx= ", D.n + D.eqidx["xz",x,z])
                append!(Var,    D.n + D.eqidx["xz",x,z])
                append!(Coeff,  -1.)
            end
        end

        # -u_{xyz} EXP* ieqs:
        for xyz = 1:D.n
            ieqn = xyz
            u_xyz = xyz
            (x,y,z) = D.xyz[xyz]
            append!(Constr, D.n + D.n + ieqn)
            println("u_rpq = ($x,$y,$z) with idx= $u_xyz")
            append!(Var, u_xyz)
            append!(Coeff, 1.)
        end

        println("co", Constr)
        println("va", Var)
        println("cof", Coeff)
        A = sparse(Constr,Var,Coeff)
        println("dimension of A", size(A))
        @show A
        return c,A,b
    end #^ make_c_b_A()


    c,A,b = make_c_A_b()
    
    var_cones = [ ( :Free, 1:n_vars) ]
    constr_cones = []

    for xyz = 1:D.n
        rqp = [ xyz, D.n + xyz, D.n + D.n + xyz ]
        append!(constr_cones, [(:ExpDual, rqp)] )
    end
    println("cstr cones ", constr_cones)
    println("var cones ", var_cones)
    
    m = MathProgBase.ConicModel(solver)
    
    MathProgBase.loadproblem!(m, c,A,b, constr_cones, var_cones)

    println("!!")
    # CPUtic()                        # TIC

    MathProgBase.optimize!(m)

    println("!!!")
    # stats_time = CPUtoq()           # TOQ


    @show MathProgBase.status(m)
    @show stats_time = MathProgBase.getsolvetime(m)


    if MathProgBase.status(m) == :Optimal || MathProgBase.status(m) == :Suboptimal || MathProgBase.status(m) == :Error
        q = Vector{Float64}(D.n)
        mu = Vector{Float64}( Mathprogbase.getdual(m) )
    
        for xyz in 1:D.n
            q[xyz] = mu[ D.n + xyz ] #or mu[3*xyz - 1] depending on the order of the reduced costs
        end

        eval_val::Float64 = InfDecomp_Base.condEntropy(D,q,big(0.))

        @show eval_val, MathProgBase.getobjval(m)
        if abs(eval_val + MathProgBase.getobjval(m)) > 1.e-3
            println("BIG DIFFERENCE IN SOLUTION VALUE!! :-( :-( :-( :-( :-( :-(\nProbably a BUG!!!!!!!!!!!!!!!!!!!!!")
        end

        stats = Solution_Stats(-eval_val,MathProgBase.getobjval(m),q,stats_time,m)
        return stats
    else
        print("Fuck!!")
        stats = Solution_Stats(999,999,[],stats_time,m)
        return stats
    end
    ; #
end #^ model_D


end #^ module
; #EOF

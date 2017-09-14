# graddesc.jl
module GradDesc

using InfDecomp_Base

using CPUTime
export my_gradient_descent

struct Solution_Stats
    optimum     :: Float64
    q           :: Vector{Float64}
    nm_pr∇      :: Float64
    max_η       :: Float64
    status      :: Symbol  # can be one of   :grad0 :step0 :iter
    time        :: Float64
end


# Compute the projection operator onto the image of the matrix B.
# Uses SVD, where singular values less than `eps` are considered 0
# Requires that B be matrix with AT LEAST AS MANY ROWS AS COLUMNS !!!!
function compute_projector( B :: AbstractMatrix{Float64}, eps::Float64=1.e-10 ) :: Matrix{Float64}
    (m,n) = size(B)
    @assert m >= n "compute_projector(B) requires that B has at least as many rows as columns"

    # SVD
    B_full = full(B)
    F = LinAlg.svdfact!(B_full)
    NZ = [ j for j in 1:n if F[:S][j] > eps ]

    # partial isometry:
    PI = view(F[:U], 1:m, NZ)

    P = PI*PI'

    @assert P*P≈P "compute_projector(B): Something went wrong -- P^2 != P :("
    @assert P'≈P  "compute_projector(B): Something went wrong -- P^T != P :("

    return P
end #^ compute_projector()


function initial_interior_point(e::My_Eval, q::Vector{Float64}) :: Void
    for x = 1:e.n_x
        p_x = 0.
        for y = 1:e.n_y
            for z = 1:e.n_z
                p_x += e.prb_xyz[x,y,z]
            end #^ for z
        end #^ for y

        for y = 1:e.n_y
            for z = 1:e.n_z
                i = e.varidx[x,y,z]
                if i>0
                    q[i] = marg_xy[x,y] * marg_xz[x,z] / p_x
                end #^ if ∃i
            end #^ for z
        end #^ for y
    end #^ for x
    ;
end

function my_gradient_descent(e::My_Eval;
                             max_iter        :: Int64    =1000000,
                             eps_grad        :: Float64  =1.e-20,
                             eps_steplength  :: Float64  =1.e-20,
                             stepfactor      :: Float64  =.1        )   :: Solution_Stats
    # Function for steplength
    steplength(q::Float64, g::Float64) = (  g > 0 ?   q/g   :  Inf  )

    CPUtic()

    # Definitions of vectors
    q_0  :: Vector{Float64} = zeros(e.n)  # standard interior feasible point (also initial q)
    q    :: Vector{Float64} = zeros(e.n)  # iterate interior feasible point
    ∇    :: Vector{Float64} = zeros(e.n)  # gradient
    pr∇  :: Vector{Float64} = zeros(e.n)  # projected gradient (onto tangent space)
    P    :: Matrix{Float64} = compute_projector(e.Gt) # projection operator (onto lin space)

    # initial solution
    initial_interior_point(e,q_0)

    # main loop
    nm_pr∇ :: Float64
    max_η  :: Float64
    status = :iter

    q .= q_0
    for iter = 1:max_iter
        # compute gradient
        InfDecomp_Base.∇f(e,∇,q,Float64(0.))

        # project gradient onto tangent space
        pr∇ .= (  P*(∇-q_0)+q_0  )

        max_η  = -1.
        nm_pr∇ = norm(pr∇)
        if nm_pr∇ ≤ eps_grad
            status = :grad0
            break
        end

        # max steplength which retains feasibility
        max_η = Inf
        for i in 1:e.n
            max_η = min(  max_η ,  steplength(q[i],pr∇[i])  )
        end

        # max eigenvalue of the Hessian
        if max_η ≤ eps_steplength
            status = :step0
            break
        end

        @show iter nm_pr∇ max_η

        q .-= (stepfactor*min(1.,max_η)) .* pr∇
    end #^ for --- main loop

    tm = CPUtoc()

    println("Terminated with")
    @show iter nm_pr∇ max_η

    return Solution_Stats(-condEntropy(e,q,Float64(0.)), q, nm_pr∇, max_η, status, tm)
end #^ my_gradient_descent()

end # module

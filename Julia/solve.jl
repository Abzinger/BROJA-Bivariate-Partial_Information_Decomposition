include("infdecomp.jl")
using InfDecomp

function my_cholesky!{FLOAT}(A::Matrix{FLOAT}) :: None
    @assert size(A)[1]==size(A)[2] "my_cholesky(): Need square matrix!"
    @assert size(A)[1] >=2         "my_cholesky(): Don't fuck with me!"

    const n::Int = size(A)[1]

    for ell in 1:n
        if A[ell,ell] <= 0
            @show ell A[ell,ell]
            @assert A[ell,ell] > 0
        end
        const d::FLOAT = sqrt(A[ell,ell])
        A[ell,ell] = d
        for k = (ell+1):n
            const c::FLOAT = A[k,ell] / d
            A[k,ell] = c
            for i = (ell+1):k
                A[k,i] -= c*A[i,ell]
            end # j
        end # k
    end # ell
    #
    for k in 1:n
        for ell in (k+1):n
            A[k,ell] = A[ell,k]
        end
    end
    return nothing;
end

import InfDecomp.My_Eval

# This function returns a pair (γ,H); the hessian is:  1/γ * H
function hessian_block{FLOAT,TFloat}(e::My_Eval, p::Vector{TFloat}, y::Int, z::Int, dummy::FLOAT) :: Tuple{FLOAT,Matrix{FLOAT}}
    @assert 1 <= y <= e.n_y
    @assert 1 <= z <= e.n_z

    P_yz  ::FLOAT = FLOAT(0.)
    for x = 1:e.n_x
        i = e.varidx[x,y,z]
        if i>0
            P_yz += p[i]
        end
    end

    H = zeros(FLOAT,e.n_x,e.n_x)

    # Start with the diagonal
    for x = 1:e.n_x
        i = e.varidx[x,y,z]
        if i>0
            P_xyz = FLOAT( p[i] )
            H[x,x] = P_yz/P_xyz - 1
        end
    end
    # Now off-diagonal.
    for x = 1:e.n_x
        for u = 1:(x-1)
            i_x = e.varidx[x,y,z]
            i_u = e.varidx[u,y,z]
            if i_x>0 && i_u>0
                counter += 1
                H[x,u] = FLOAT(-1)
                H[u,x] = FLOAT(-1)
            end
        end
    end


    return (P_yz,H)
    ;
end


The Quadratic Program
---------------------

[ a ]     [ H   -A^T ] [ u ]
[   ] ==  [          ] [   ]
[ b ]     [ A     0  ] [ v ]

For Schur-complement method, need to solve:

  A U (U^T H U)^{-1} U^T A^T v == c

where U is the isometry taking im(A^T) into the full space; H is invertible on that space.

The equation

   AU w = d

has a unique solution w in im(A^T).

Then solve A^T v = H w    for v.

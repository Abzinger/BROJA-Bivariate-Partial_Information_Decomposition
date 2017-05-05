using InfDecomp

function my_cholesky!{FLOAT}(A::Matrix{FLOAT})
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
    ;
end


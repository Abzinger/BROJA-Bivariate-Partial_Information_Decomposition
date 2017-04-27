using Cuba

denom(Σinv::Matrix{Float64}) :: BigFloat  =  sqrt( (2*π)^3 / det(Σinv) )
function gaussint(Σinv::Matrix{Float64}, x::Vector{Float64}) :: BigFloat
    t = - BigFloat( (x'*Σinv*x)[1] )/2
    return exp(t)
    ;
end

function integrate(Σinv::Matrix{Float64}, a::Vector{Float64}, b::Vector{Float64}, alg, n::Int=100) :: Float64
    function fun(x::Vector{Float64}, f::Vector{Float64})
        tmp = BigFloat(0)
        y = a + x.*d
        for di in -1:+1
            for dj in -1:+1
                for dk in -1:+1
                    v = gaussint(Σinv, y+[di/1000000,dj/1000000,dk/1000000].*d)
                    di!=0 || dj!=0 || dk!=0 || ( v*=4 )
                    tmp += v
                end
            end
        end
        f[1] = Float64( tmp )
    end

    divi = denom(Σinv)
    d = b-a

    if alg==:vegas
        integral,error,probability,neval,fail,nregions = vegas(fun,3,1)
        integral[1] /= d[1]*d[2]*d[3]*divi*12
        error[1]    /= d[1]*d[2]*d[3]*divi*12
    elseif alg==:suave
        integral,error,probability,neval,fail,nregions = suave(fun,3,1)
        integral[1] /= d[1]*d[2]*d[3]*divi*12
        error[1]    /= d[1]*d[2]*d[3]*divi*12
    elseif alg==:divonne
        integral,error,probability,neval,fail,nregions = divonne(fun,3,1)
        integral[1] /= d[1]*d[2]*d[3]*divi*12
        error[1]    /= d[1]*d[2]*d[3]*divi*12
    elseif alg==:cuhre
        integral,error,probability,neval,fail,nregions = cuhre(fun,3,1)
        integral[1] /= d[1]*d[2]*d[3]*divi*12
        error[1]    /= d[1]*d[2]*d[3]*divi*12
    elseif alg==:grid
        x = Vector{Float64}([0.,0.,0.])
        avg = BigFloat(0)
        for i = 0:(n-1)
            x[1] = ( i*b[1] + (n-i)*a[1] )/n
            for j = 0:(n-1)
                x[2] = ( j*b[2] + (n-j)*a[2] )/n
                for k = 0:(n-1)
                    x[3] = ( k*b[3] + (n-k)*a[3] )/n
                    p = gaussint(Σinv, x) / d
                    avg +=  p / n^3
                end
            end
        end
        integral,error,probability,neval,fail,nregions = avg,0,0,n,0,0
    else
        @assert false "Wrong Algorithm symbol"
    end
    @show integral error probability neval fail nregions
    return integral[1]
end


function make_apxgaussian_pdf(Σ::Matrix{Float64}, blocklength::Float64, zero::Float64=1.e-50, alg=:divonne, n::Int = 100) :: Dict{Tuple{Float64,Float64,Float64},Float64}
    @assert n>=4 "make_apxgaussian_pdf(): n>=4 needed"
    Σinv = inv(Σ)

    pdf = Dict{Tuple{Float64,Float64,Float64},Float64}()

    prc = precision(BigFloat)
    setprecision(BigFloat,512)

    boundary = Array{Tuple{Int,Int,Int}}(1)
    empty!(boundary)
    # Do central box, initialize boundary.
    i_,j_,k_ = (0,0,0)
    a = [ i_   *blocklength,   j_   *blocklength,   k_   *blocklength]
    b = [(i_+1)*blocklength,  (j_+1)*blocklength,  (k_+1)*blocklength]
    c = (a.+b)./2
    p = Float64( integrate(Σinv,a,b,alg,n) )
    if p > zero
        pdf[ (c[1],c[2],c[3]) ] = p
        push!(boundary, (i_,j_,k_) )
        # @show (i_,j_,k_) p
    end

    while !isempty(boundary)
        i,j,k = pop!(boundary)
        # @show (i,j,k) length(boundary)
        for di in -1:+1
            for dj in -1:+1
                for dk in -1:+1
                    i_,j_,k_ = (i+di,j+dj,k+dk)
                    a = [ i_   *blocklength,   j_   *blocklength,   k_   *blocklength]
                    b = [(i_+1)*blocklength,  (j_+1)*blocklength,  (k_+1)*blocklength]
                    c = (a.+b)./2
                    if (c[1],c[2],c[3]) ∉ keys(pdf)
                        p = Float64( integrate(Σinv,a,b,alg,n) )
                        # @show c p
                        if p > zero
                            pdf[ (c[1],c[2],c[3]) ] = p
                            push!(boundary, (i_,j_,k_) )
                            # @show (i_,j_,k_) p
                        end
                    end
                end
            end
        end
    end

    setprecision(BigFloat,prc)
    return pdf
    ;
end

random_Σ() :: Matrix{Float64} = begin A=ones(3,3)-2*rand(3,3); return A'*A ; end

include("infdecomp.jl")
using InfDecomp

list_of_functions = [
                     :rdn,
                     :unq,
                     :xor,
                     :and,
                     :rdnxor,
                     :rdnunqxor,
                     :xorand,
#                     :discr_normal
                     ]


function sample(fun::Symbol, Pr_noise::Float64, sz)
    if fun==:rdn # X=Y=Z
        @assert sz ≥ 2  "Need size ≥ 2 for rdn"
        y=rand(0:sz-1)
        z=y
        x=y
        if rand() ≤ Pr_noise
            x = mod(x + rand([+1,-1])  ,sz)
        end
        return (x,y,z)
    elseif fun==:unq # (Y,Z),Y,Z
        @assert sz ≥ 2  "Need size ≥ 2 for unq"
        y=rand(0:sz-1)
        z=rand(0:sz-1)
        x=(y,z)
        if rand() ≤ Pr_noise
            while x==(y,z)
                x = (  mod(y+rand([+1,0,-1]),sz), mod(z+rand([+1,0,-1]),sz)  )
            end
        end
        return ( x,y,z )
    elseif fun==:xor # (Y xor Z),Y,Z
        @assert ispow2(sz)  "Need size to be a power of 2 for xor"
        y=rand(0:sz-1)
        z=rand(0:sz-1)
        x= y $ z
        if rand() ≤ Pr_noise
            l = Int(round(log2(sz)))
            b = rand(0:l-1)
            N = 2^b
            x $= N
        end
        return ( x,y,z )
    elseif fun==:and # (Y and Z),Y,Z
        @assert ispow2(sz)  "Need size to be a power of 2 for and"
        y=rand(0:sz-1)
        z=rand(0:sz-1)
        x= y & z
        if rand() ≤ Pr_noise
            l = Int(round(log2(sz)))
            b = rand(0:l-1)
            N = 2^b
            x $= N
        end
        return ( x,y,z )
    elseif fun==:rdnxor # (  (U xor V,W),  (U,W), (V,W) )
        @assert ispow2(sz)  "Need size to be a power of 2 for rdnxor"
        u=rand(0:sz-1)
        v=rand(0:sz-1)
        w=rand(0:sz-1)
        y = (u,w)
        z = (u,w)
        x = ( u $ v , w)
        if rand() ≤ Pr_noise
            l = Int(round(log2(sz)))
            b = rand(0:l-1)
            N = 2^b
            r = rand()
            if r ≤ 2/5
                x = ( u $ v $ N, w)
            elseif r ≤ 4/5
                x = ( u $ v, w $ N)
            else
                b = rand(0:l-1)
                M = 2^b
                x = ( u $ v $ M, w $ N)
            end
        end
        return ( x,y,z )
    elseif fun==:xorand
        @assert ispow2(sz)  "Need size to be a power of 2 for and"
        y=rand(0:sz-1)
        z=rand(0:sz-1)
        x= (y$z, y&z)
        if rand() ≤ Pr_noise
            l = Int(round(log2(sz)))
            b = rand(0:l-1)
            N = 2^b
            if rand() ≤ 1/2
                x = (y$z $ N, y&z)
            else
                x = (y$z, y&z $ N)
            end
        end
        return ( x,y,z )
    elseif fun==:rdnunqxor # (
        @assert typeof(sz) <: Tuple  "Need size to be tuple for rdnunqxor"
        @assert length(sz) == 2      "Need size to be pair for rdnunqxor"
        @assert ispow2(sz[1])        "Need sizes to be powers of 2 for rdnunqxor"
        @assert sz[1]>=4             "Need 1st size to be ≥ 2 for rdnunqxor"
        l  = Int(round(log2(sz[1])))
        l₁ = Int(floor(l/2))
        l₂ = Int(ceil(l/2))
        s₁ = 2^l₁
        s₂ = 2^l₂
        u₁ = rand(0:s₁-1)
        u₂ = rand(0:s₂-1)
        v₁ = rand(0:s₁-1)
        v₂ = rand(0:s₂-1)
        w  = rand(0:sz[2]-1)
        x = (u₁$v₁,u₂,v₂,w)
        y = (u₁,u₂,w)
        z = (v₁,v₂,w)
        if rand() ≤ Pr_noise
            r = rand()
            if r ≤ 1/5
                u₁ = rand(0:s₁-1)
            elseif r ≤ 2/5
                u₂ = rand(0:s₂-1)
            elseif r ≤ 3/5
                v₁ = rand(0:s₁-1)
            elseif r ≤ 4/5
                v₂ = rand(0:s₂-1)
            else
                w  = rand(0:sz[2]-1)
            end
            x = (u₁$v₁,u₂,v₂,w)
        end
        return ( x,y,z )
    # elseif fun==:discr_normal
    #     blah
    else
        error("sample(): Unknown function, $fun")
    end
    ;
end


function make_PDF(fun::Symbol, n_samples::Int, Pr_noise::Float64, sz, FLOAT::DataType=Float64)
    t = sample(fun,0.,sz)
    T = typeof(t)
    @assert n_samples >= 1
    @assert T <: Tuple "Something went wrong, I didn't get a tuple!"
    @assert length(t) == 3 "Something went wrong, I didn't get a 3-tuple!"

    pdf =  Dict{T,FLOAT}()

    for iter = 1:n_samples
        t = sample(fun,Pr_noise,sz)
        pdf[t] = get!(pdf,t,FLOAT(0)) + FLOAT(1)
    end

    for t in keys(pdf)
        pdf[t] /= n_samples
    end

    return pdf
    ;
end # make_PDF()

function test()
    for f in list_of_functions
        if f ≠ :rdnunqxor
            q = make_PDF(f,100,.05,16)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,8)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,4)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,2)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,16)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,8)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,4)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,2)
            d,e = InfDecomp.create_stuff(q)
        else
            q = make_PDF(f,100,.05,(16,3))
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,(8,3))
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,(4,3))
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,(16,5))
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,(8,5))
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,(4,5))
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,(16,3))
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,(8,3))
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,(4,3))
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,(16,5))
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,(8,5))
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,(4,5))
            d,e = InfDecomp.create_stuff(q)
        end
        if f ≠ :rdnunqxor
            q = make_PDF(f,100,.05,16,BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,8,BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,4,BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,2,BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,16,BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,8,BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,4,BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,2,BigFloat)
            d,e = InfDecomp.create_stuff(q)
        else
            q = make_PDF(f,100,.05,(16,3),BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,(8,3),BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,(4,3),BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,(16,5),BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,(8,5),BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.05,(4,5),BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,(16,3),BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,(8,3),BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,(4,3),BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,(16,5),BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,(8,5),BigFloat)
            d,e = InfDecomp.create_stuff(q)
            q = make_PDF(f,100,.95,(4,5),BigFloat)
            d,e = InfDecomp.create_stuff(q)
        end
    end
end


; # EOF

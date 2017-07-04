include("gaussian.jl")

list = [1.0, 0.75, 0.5, 0.4, 0.25]
for b in 16:20
    Σ = random_Σ()
    writedlm("Data/gauss-$b.matrix", Σ)
    for s in list
        pdf = make_apxgaussian_pdf(Σ,s)
        open("Data/gauss-$b-$s.dens", "a") do dfile
            print(dfile,"[")
            i = 1
            for (k, v) in pdf
                j = 1
                if i < length(pdf)
                    print(dfile,"[[")
                    for w in k
                        if j < length(k)
                            print(dfile,"$w, ")
                        else
                            print(dfile,"$w")
                        end
                        j +=1
                    end
                    print(dfile,"], $v], ")
                else
                    print(dfile,"[[")
                    for w in k
                        if j < length(k) 
                            print(dfile,"$w, ")
                        else
                            print(dfile,"$w")
                        end
                        j += 1
                    end
                    print(dfile,"], $v]")
                end
                i += 1
            end
            print(dfile,"]")
        end
    end
end


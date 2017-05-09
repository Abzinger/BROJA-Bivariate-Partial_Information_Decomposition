include("infdecomp.jl")
include("read_pdf.jl")
using InfDecomp
using Read_PDF
using MathProgBase 
using Mosek
using Ipopt
using NLopt
using KNITRO
#-------------------------------
#READS P (CONTROL) PROBABILITY
#-------------------------------
function read_p(w::String)
	 pdf = Dict()
	 #pdf = read_pdf(w)
	 pdf = read_pdf("$w")
	 return pdf
	;
end

#-------------------------------------------------
#COMPUTES TOTAL VARIATION FROM TRUE DISTRIBUTION
#-------------------------------------------------
function tot_var{T1,T2,T3}(pdf::Dict{Tuple{T1,T2,T3},Float64}, true_pdf::Dict{Tuple{T1,T2,T3},Float64})
    a  = Float64[]
    b  = Float64[]
    u = 0
    w = Float64[]
    for (k,v) in pdf
        u = abs(v - get(true_pdf,k,0))
        push!(a, u)
    end
    for (k,v) in true_pdf
        if get(pdf,k,0) == 0
            push!(b, v)
        end
    end
    w = append!(a,b)
    return 0.5 * norm(w, 1)
end
#-------
#RUN IT
#-------
function test(solver::Symbol,w1::String, w2::String)
    p = read_p(w1)
    true_p = read_p(w2)
    if (solver == :Mosek)
        println("Start optimization with Mosek...")
	sd,myeval,model = InfDecomp.do_it(p, Mosek.MosekSolver(MSK_IPAR_INTPNT_MULTI_THREAD=0, MSK_IPAR_INTPNT_MAX_ITERATIONS=500))
    elseif (solver == :Ipopt)
        println("Start optimization with Ipopt...")
	sd,myeval,model = InfDecomp.do_it(p,Ipopt.IpoptSolver())
    elseif (solver == :NLopt_LD_SLSQP)
        println("Start optimization with NLopt/LD_SLSQP...")
	sd,myeval,model = InfDecomp.do_it(p,NLopt.NLoptSolver(algorithm=:LD_SLSQP))
    elseif (solver == :Knitro)
        println("Start optimization with Knitro...")
	sd,myeval,model = InfDecomp.do_it(p,KNITRO.KnitroSolver())
    end
    feasstats = InfDecomp.check_feasibility(model,myeval)
    s = tot_var(p,true_p)
    open("feas_stats.csv", "a") do ffile
        print(ffile,"# filename")
        for i in fieldnames(feasstats)
            print(ffile,",\t",i)
        end
        print(ffile,", Tot_var\n",w1)
        for i in fieldnames(feasstats)
            print(ffile,",\t",getfield(feasstats,i))
        end
        print(ffile,", $s\n")
    end
    return sd,myeval,model
end


; # EOF

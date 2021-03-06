include("infdecomp.jl")
include("read_pdf.jl")
using InfDecomp
using Read_PDF
using MathProgBase
using Mosek
using Ipopt
using NLopt
using KNITRO
using PyCall
unshift!(PyVector(pyimport("sys")["path"]), "")
@pyimport cvxopt_solve
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

#-------------------------------------------------
# Solutions & Statistics for Cvxopt
#-------------------------------------------------
# type Solution_and_Stats
#     var_num             :: Int64
#     x_sz                :: Int64
#     y_sz                :: Int64
#     z_sz                :: Int64
#     status              :: Symbol
#     obj_val             :: BigFloat
#     q_nonneg_viol       :: BigFloat
#     marginals_1         :: BigFloat
#     marginals_2         :: BigFloat
#     marginals_Inf       :: BigFloat
#     mu_nonneg_viol      :: BigFloat
#     complementarity_max :: BigFloat
#     complementarity_sum :: BigFloat
#     # MI_X_YZ           :: BigFloat
#     CI                  :: BigFloat
#     SI                  :: BigFloat
#     UI_Y                :: BigFloat
#     UI_Z                :: BigFloat
#     opt_time            :: String
#     # entropy_X           :: BigFloat
# end

Solution_and_Stats = InfDecomp.Solution_and_Stats

#-------
#RUN IT
#-------
function test(solver::Symbol,w1::String, w2::String, tmpFloatDatatype::DataType)
    p = read_p(w1)
    true_p = read_p(w2)
    time_l = 10000.0
    if (solver != :Cvxopt)
        if (solver == :Mosek)
            InfDecomp.set_copy_sol_behaviour(true)
            println("Start optimization with Mosek...")
	    sd,myeval,model = InfDecomp.do_it(p, Mosek.MosekSolver(MSK_IPAR_INTPNT_MULTI_THREAD=0, MSK_IPAR_INTPNT_MAX_ITERATIONS=500, MSK_DPAR_OPTIMIZER_MAX_TIME=time_l),tmpFloatDatatype)
        elseif (solver == :Ipopt)
            InfDecomp.set_copy_sol_behaviour(false)
            println("Start optimization with Ipopt...")
	    sd,myeval,model = InfDecomp.do_it(p,Ipopt.IpoptSolver(max_cpu_time=time_l),tmpFloatDatatype)
            # elseif (solver == :NLopt_LD_SLSQP)
            #     println("Start optimization with NLopt/LD_SLSQP...")
            #     sd,myeval,model = InfDecomp.do_it(p,NLopt.NLoptSolver(algorithm=:LD_SLSQP))
        elseif (solver == :Knitro_Ip)
            InfDecomp.set_copy_sol_behaviour(false)
            println("Start optimization with Knitro/Interior_pt...")
	    sd,myeval,model = InfDecomp.do_it(p,KNITRO.KnitroSolver(algorithm=1,maxtime_cpu=time_l),tmpFloatDatatype)
        elseif (solver == :Knitro_IpCG)
            InfDecomp.set_copy_sol_behaviour(false)
            println("Start optimization with Knitro/Interior_pt_CG...")
	    sd,myeval,model = InfDecomp.do_it(p,KNITRO.KnitroSolver(algorithm=2,maxtime_cpu=time_l),tmpFloatDatatype)
        elseif (solver == :Knitro_AS)
            InfDecomp.set_copy_sol_behaviour(false)
            println("Start optimization with Knitro/Active_Set...")
	    sd,myeval,model = InfDecomp.do_it(p,KNITRO.KnitroSolver(algorithm=3,maxtime_cpu=time_l),tmpFloatDatatype)
        elseif (solver == :Knitro_SQP)
            InfDecomp.set_copy_sol_behaviour(false)
            println("Start optimization with Knitro/SQP...")
	    sd,myeval,model = InfDecomp.do_it(p,KNITRO.KnitroSolver(algorithm=4,maxtime_cpu=time_l),tmpFloatDatatype)
        end
        feasstats = InfDecomp.check_feasibility(model,myeval,solver)
        s = tot_var(p,true_p)
        open("feas_stats_64.csv", "a") do ffile
            print(ffile,"# filename, Solver")
            for i in fieldnames(feasstats)
                print(ffile,",\t",i)
            end
            if (solver ==:Mosek)
                print(ffile,", Tot_var\n",w1,",Mosek")
            elseif (solver ==:Ipopt)
                print(ffile,", Tot_var\n",w1,",Ipopt")
            elseif (solver ==:Knitro_Ip)
                print(ffile,", Tot_var\n",w1,",Knitro_Ip")
            elseif (solver ==:Knitro_IpCG)
                print(ffile,", Tot_var\n",w1,",Knitro_IpCG")
            elseif (solver ==:Knitro_AS)
                print(ffile,", Tot_var\n",w1,",Knitro_AS")
            elseif (solver ==:Knitro_SQP)
                print(ffile,", Tot_var\n",w1,",Knitro_SQP")
            end
            for i in fieldnames(feasstats)
                print(ffile,",\t",getfield(feasstats,i))
            end
            print(ffile,", $s\n")
        end
        return sd,myeval,model
    else
        println("Start optimization with Cvxopt...")
        pdf = PyDict(p)
        fstat  = Solution_and_Stats( 0,0,0,0," ",  0,0,0,0,0,0,0,0,0,0,0,0," ")
        fstat.var_num, fstat.x_sz, fstat.y_sz, fstat.z_sz, fstat.status, fstat.obj_val, fstat.q_nonneg_viol, fstat.marginals_1, fstat.marginals_2, fstat.marginals_Inf, fstat.mu_nonneg_viol, fstat.complementarity_max, fstat.complementarity_sum, fstat.CI, fstat.SI, fstat.UI_Y, fstat.UI_Z = cvxopt_solve.solve_PDF(p,time_l*1000)
        fstat.opt_time = -1

        s = tot_var(p,true_p)
        open("feas_stats_64_gu_9.csv", "a") do ffile
            print(ffile,"# filename, Solver")
            for i in fieldnames(fstat)
                print(ffile,",\t",i)
            end
            print(ffile,", Tot_var\n",w1,",Cvxopt")
            for i in fieldnames(fstat)
                print(ffile,",\t",getfield(fstat,i))
            end
            print(ffile,", $s\n")
        end
        return 0
    end
end


; # EOF

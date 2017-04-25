include("infdecomp.jl")
include("read_pdf.jl")
using InfDecomp
using Read_PDF
using MathProgBase 
using Mosek
using Ipopt

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
#----------------------
#TABULATE THE RESULTS
#----------------------
function create_csv(comp_1::BigFloat, comp_inf::BigFloat, dist_1::Float64, dist_2::Float64, dist_inf::Float64)
    open("table.csv", "a") do f
        println(f,"$comp_1,$comp_inf,$dist_1,$dist_2,$dist_inf")
    end
    ;
end

#-------
#RUN IT
#-------
function test(solver::Symbol,w::String)
    p = read_p(w)
    if (solver == :Mosek)
        println("Start optimization with Mosek...")
	sd,myeval,model = InfDecomp.do_it(p, Mosek.MosekSolver())
        comp_1, comp_inf = InfDecomp.KKT_feas(model, myeval)
        println("Complementary slackness computed...",)
        dist_1, dist_2, dist_inf = InfDecomp.marg_dist(model, myeval)
        println("Primal Feasibility violations computed...")
        create_csv(comp_1,comp_inf,dist_1,dist_2,dist_inf)
        println("Table is created...")
    elseif (solver == :Ipopt)
        println("Start optimization with Ipopt...")
	sd,myeval,model = InfDecomp.do_it(p,Ipopt.IpoptSolver())
        comp_1, comp_inf = InfDecomp.KKT_feas(model, myeval)
        println("Complementary slackness computed...",)
        dist_1, dist_2, dist_inf = InfDecomp.marg_dist(model, myeval)
        println("Primal Feasibility violations computed...")
        create_csv(comp_1,comp_inf,dist_1,dist_2,dist_inf)
        println("Table is created...")
    end

end


; # EOF

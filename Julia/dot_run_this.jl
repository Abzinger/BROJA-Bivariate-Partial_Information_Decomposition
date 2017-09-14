include("infdecomp_base.jl")
include("expcone.jl")
include("graddesc.jl")
include("infdecomp.jl")
include("read_pdf.jl")
include("dot_try.jl")
using InfDecomp
# solvers_list = [:Mosek, :Ipopt, :Knitro_Ip, :Knitro_IpCG, :Knitro_AS, :Knitro_SQP, :Cvxopt]

#solvers_list = [:Mosek]

tmpFloatDatatype = Float64


if length(ARGS) != 2
    print(
    """Use w/ two arguments: solver & instance.
    Solver:  one of ECOS_L, Mosek, GD.
    instance: the name of a file in ../PDFs/Data/
    """)
    exit(1)
end

@show ARGS

solver = nothing

if     uppercase(ARGS[1]) == "ECOS_L"
    solver = :ECOS_L
elseif uppercase(ARGS[1]) == "MOSEK"
    solver = :Mosek
elseif uppercase(ARGS[1]) == "GD"
    solver = :My_GradDesc
else
    println("I don't recognize the solver (only ECOS_L, Mosek, GD are allowed)")
    exit(1)
end

instance = ARGS[2]

test(solver, "AND-0.1-10000.dens", "AND-0.1-10000.dens",  tmpFloatDatatype)

stats_file = "dot_stats.csv"

# open(stats_file, "a") do ffile
#     for i in fieldnames(InfDecomp.Solution_and_Stats)
#         print(ffile,"\t",i,", ")
#     end
#     println(ffile,"")
# end

println("########################################################################################################################################################################################################")
println("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
println("Starting computations of the distribution: ", instance," with ", solver)
println("")
test(solver, instance, instance,  tmpFloatDatatype, stats_file=stats_file)
println("")
println("Finishing computations of the distribution: ", instance," with ", solver)
println("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

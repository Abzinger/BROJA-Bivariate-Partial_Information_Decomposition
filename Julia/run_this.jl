include("try.jl")


# solvers_list = [:Mosek, :Ipopt, :Knitro_Ip, :Knitro_IpCG, :Knitro_AS, :Knitro_SQP, :Cvxopt]

solvers_list = [:Mosek]

list= String[]
for i in 1:5
    push!(list, "Data/gauss-$i")
end
size =[1,0.75,0.5,0.4]
list_instances=String[]
for l in 1:length(list)
    for s in size
        push!(list_instances, string(list[l],"-$s.dens"))
    end#^n
end#^s

for solver in solvers_list
    for instance in list_instances
        println("///////////////////////////")
        println("Starting computations of the distribution: ", instance," with ", solver)
        println("///////////////////////////")
        test(solver, instance, instance)
        println("///////////////////////////")
        println("Finishing computations of the distribution: ", instance," with ", solver)
        println("///////////////////////////")
    end
end


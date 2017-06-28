include("try.jl")

 true_list = ["../PDFs/Data/XOR-0.0-0.dens","../PDFs/Data/AND-0.0-0.dens","../PDFs/Data/UNQ-0.0-0.dens","../PDFs/Data/RDN-0.0-0.dens","../PDFs/Data/XORAND-0.0-0.dens","../PDFs/Data/RDNXOR-0.0-0.dens","../PDFs/Data/RDNUNQXOR-0.0-0.dens"]
# true_list = ["../PDFs/Data/XOR-0.0-0.dens"]
# true_list = ["../PDFs/Data/AND-0.0-0.dens"]
# true_list = ["../PDFs/Data/UNQ-0.0-0.dens"]
# true_list = ["../PDFs/Data/RDN-0.0-0.dens"]
# true_list = ["../PDFs/Data/XORAND-0.0-0.dens"]
# true_list = ["../PDFs/Data/RDNXOR-0.0-0.dens"]
# true_list = ["../PDFs/Data/RDNUNQXOR-0.0-0.dens"]

 list = ["../PDFs/Data/XOR-0.05-100000.dens","../PDFs/Data/AND-0.05-100000.dens","../PDFs/Data/UNQ-0.05-100000.dens","../PDFs/Data/RDN-0.05-100000.dens","../PDFs/Data/XORAND-0.05-100000.dens","../PDFs/Data/RDNXOR-0.05-100000.dens","../PDFs/Data/RDNUNQXOR-0.05-100000.dens"]
# list = ["../PDFs/Data/XOR-0.1-10000.dens"]
# list = ["../PDFs/Data/AND-0.1-10000.dens"]
# list = ["../PDFs/Data/UNQ-0.1-10000.dens"]
# list = ["../PDFs/Data/RDN-0.1-10000.dens"]
# list = ["../PDFs/Data/XORAND-0.1-10000.dens"]
# list = ["../PDFs/Data/RDNXOR-0.1-10000.dens"]
# list = ["../PDFs/Data/RDNUNQXOR-0.1-10000.dens"]
list_new=String[]
for s in 1:length(list)
    for n in 1:101
        if n == 1
            push!(list_new, list[s])
        elseif n == 101
            continue
        end
        r = n - 1
        push!(list_new, string(list[s],"-$r"))
    end#^n
end#^s

# for s in 1:length(list)
#     for n in 1:10
#         push!(list_new, string(list[s],"-$n"))
#     end#^n
# end#^s

# solvers_list = [:Mosek, :Ipopt, :Knitro_Ip, :Knitro_IpCG, :Knitro_AS, :Knitro_SQP, :Cvxopt]
solvers_list = [:Ipopt]
tmpFloatDatatype = Float64
for s in 1:length(list_new)
    for solver in solvers_list
        println("///////////////////////////")
        println("Starting computations of the distribution: ", list_new[s], " with ", solver, " and prc ", tmpFloatDatatype)
        println("///////////////////////////")
        test(solver,list_new[s], true_list[convert(Int32,ceil(s/101))], tmpFloatDatatype)
        println("///////////////////////////")
        println("Finishing computations of the distribution: ", list_new[s]," with ", solver, " and prc ", tmpFloatDatatype)
        println("///////////////////////////")
    end
end

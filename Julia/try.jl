include("infdecomp.jl")
include("read_pdf.jl")
using InfDecomp
using Read_PDF


function try_dist(w::String)
	 pdf = Dict()
	 #pdf = read_pdf(w)
	 pdf = read_pdf("Data/$w")
	 return pdf	
	;
end
function test(i::Int64,w::String)
	 q = try_dist(w)
if (i == 0)
	InfDecomp.do_it(q,Mosek.MosekSolver())
elseif (i == 1)
	InfDecomp.do_it(q,Ipopt.IpoptSolver())
end
	
end



; # EOF
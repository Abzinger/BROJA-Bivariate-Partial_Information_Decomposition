include("infdecomp.jl")
include("read_pdf.jl")
using InfDecomp
using Read_PDF

using Mosek
using Ipopt


function try_dist(w::String)
	 pdf = Dict()
	 #pdf = read_pdf(w)
	 pdf = read_pdf("$w")
	 return pdf
	;
end
function test(solver::Symbol,w::String)
    q = try_dist(w)
    if (solver == :Mosek)
	InfDecomp.do_it(q,Mosek.MosekSolver())
    elseif (solver == :Ipopt)
	InfDecomp.do_it(q,Ipopt.IpoptSolver())
    end

end



; # EOF

#read_pdf.jl
module Read_PDF
export read_pdf

using JSON

function to_tuple(X) :: Any
    if typeof(X) <: Array
        tmp = [ to_tuple(x) for x in X ]
        x = (tmp[1:end]...)
    else
        return X
    end
end


function array_to_dict(a) :: Dict{Tuple{Any,Any,Any}, Float64}
    local d = Dict{Tuple{Any,Any,Any}, Float64}()
    for p in a
        t = Tuple{Any,Any,Any}( to_tuple( p[1] ) )
        d[t] = Float64( p[2] )
    end
    return d
end


function read_pdf(filename::String) ::  Dict{Tuple{Any,Any,Any}, Float64}
    a = JSON.parsefile(filename)
    return array_to_dict(a)
end


end
